#pragma once

#include <Cosserat/engine/CosseratIntrinsicState.h>
#include <liegroups/SO3.h>

#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/SingleLink.h>
#include <sofa/type/Mat.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

namespace sofa::component::cosserat::engine {

/**
 * @brief Force field for Cosserat beams using the staggered (Painless) discretisation.
 *
 * Implements elastic forces and the tangent stiffness for the staggered Cosserat
 * beam described in:
 *
 *   Romanyà-Serrasolsas, Casafranca, Otaduy —
 *   "Painless Differentiable Rotation Dynamics", ACM ToG / SIGGRAPH 2025.
 *   DOI: 10.1145/3730944
 *
 * ## DOF layout (stored in a linked CosseratIntrinsicState)
 *
 *   x_0, x_1, …, x_N        ∈ R³   (N+1 node positions)
 *   R_0, R_1, …, R_{N-1}    ∈ SO3  (N  segment orientations, at edge midpoints)
 *
 * ## Elastic energies
 *
 * ### Linear (stretch + shear) — one contribution per segment i
 * @code
 *   Γ_i   = R_i⁻¹ · (x_{i+1} − x_i) / h_i − e₁
 *   E_lin = Σ_i  h_i/2 · Γ_i^T · K_L · Γ_i
 *   K_L   = diag(EA, GA, GA)
 * @endcode
 *
 * ### Angular (bending + torsion) — one contribution per interior node i
 * @code
 *   φ_i   = log(R_{i-1}^T · R_i)
 *   h̃_i   = (h_{i-1} + h_i) / 2          (dual-edge length)
 *   Ω_i   = φ_i / h̃_i
 *   E_ang = Σ_i  h̃_i/2 · Ω_i^T · K_A · Ω_i
 *   K_A   = diag(GJ, EIy, EIz)
 * @endcode
 *
 * ## Forces (negative gradient of energy)
 *
 * ### On position DOFs (world frame)
 * @code
 *   f(x_i)     += +R_i · K_L · Γ_i
 *   f(x_{i+1}) += −R_i · K_L · Γ_i
 * @endcode
 *
 * ### On SO3 DOFs — torques in the body frame of each segment
 * Using the identity J_r^{-T}(φ) = J_r^{-1}(−φ):
 * @code
 *   τ(R_i)     += −J_r^{-1}(−φ_i) · K_A · Ω_i      (from node i)
 *   τ(R_{i-1}) += +J_r^{-1}( φ_i) · K_A · Ω_i      (from node i)
 * @endcode
 * accumulated over all interior nodes i ∈ {1, …, N−1}.
 *
 * @note addDForce and addKToMatrix implement the **linear** stiffness blocks.
 *       The angular (geometric + material) stiffness blocks are marked TODO
 *       and will be added once the linear path is validated.
 */
class PainlessBeamForceField : public sofa::core::behavior::BaseForceField {
   public:
    SOFA_CLASS(PainlessBeamForceField, sofa::core::behavior::BaseForceField);

    using Vec3d    = sofa::type::Vec3d;
    using Mat3x3d  = sofa::type::Mat3x3d;
    using SO3      = sofa::component::cosserat::liegroups::SO3<double>;
    using VecVec3d = sofa::type::vector<Vec3d>;

    // ── Stiffness parameters ───────────────────────────────────────────────────
    sofa::core::objectmodel::Data<double> d_EA;   ///< Axial stiffness  E·A
    sofa::core::objectmodel::Data<double> d_GA;   ///< Shear stiffness  G·A  (both transverse directions)
    sofa::core::objectmodel::Data<double> d_GJ;   ///< Torsion stiffness G·J
    sofa::core::objectmodel::Data<double> d_EIy;  ///< Bending stiffness E·Iy (y-axis)
    sofa::core::objectmodel::Data<double> d_EIz;  ///< Bending stiffness E·Iz (z-axis)

    // ── Link to the mechanical state ───────────────────────────────────────────
    sofa::core::objectmodel::SingleLink<CosseratIntrinsicState> l_state;

    // ── Read-only outputs (debugging / visualisation / post-processing) ────────
    /// Elastic forces on N+1 position DOFs, expressed in world frame.
    sofa::core::objectmodel::Data<VecVec3d> d_nodalForces;
    /// Elastic torques on N SO3 DOFs, expressed in the body frame of each segment.
    sofa::core::objectmodel::Data<VecVec3d> d_segmentTorques;

    // ─────────────────────────────────────────────────────────────────────────
    PainlessBeamForceField();
    ~PainlessBeamForceField() override = default;

    void init() override;
    void reinit() override;

    // ── BaseForceField interface ───────────────────────────────────────────────

    /**
     * @brief Compute elastic forces and accumulate them.
     *
     * Results are written to d_nodalForces (world frame) and
     * d_segmentTorques (body frame).  Full integration into the SOFA
     * solver requires CosseratIntrinsicState to expose standard VecDeriv
     * accessors; that step is tracked as a separate TODO.
     */
    void addForce(const sofa::core::MechanicalParams* mparams,
                  sofa::core::MultiVecDerivId f,
                  sofa::core::ConstMultiVecCoordId x,
                  sofa::core::ConstMultiVecDerivId v) override;

    /**
     * @brief Compute differential elastic forces  df = −kFactor · K_lin · dx.
     *
     * Implements the **linear** stiffness contribution (stretch + shear).
     * The angular (bending/torsion) stiffness differential is TODO.
     */
    void addDForce(const sofa::core::MechanicalParams* mparams,
                   sofa::core::MultiVecDerivId df,
                   sofa::core::ConstMultiVecDerivId dx) override;

    /**
     * @brief Assemble the tangent stiffness into the global matrix.
     *
     * Adds the 3×3 linear stiffness blocks:
     *   K_{ii}, K_{i+1,i+1}  += +R_i · K_L · R_i^T / h_i
     *   K_{i,i+1}            += −R_i · K_L · R_i^T / h_i   (and transpose)
     *
     * The angular stiffness blocks are TODO.
     */
    void addKToMatrix(const sofa::core::MechanicalParams* mparams,
                      const sofa::core::behavior::MultiMatrixAccessor* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams* mparams,
                             sofa::core::ConstMultiVecCoordId x) const override;

   private:
    // ── Internal helpers ───────────────────────────────────────────────────────

    /// Returns K_L = diag(EA, GA, GA) in the body frame.
    Mat3x3d buildK_L() const;

    /// Returns K_A = diag(GJ, EIy, EIz) in the body frame.
    Mat3x3d buildK_A() const;

    /**
     * @brief Core force/torque computation from current state.
     *
     * @param[out] f_nodes    Forces on N+1 position DOFs (world frame).
     * @param[out] tau_segs   Torques on N SO3 DOFs (body frame).
     * @return Elastic potential energy.
     */
    double computeForcesAndTorques(VecVec3d& f_nodes, VecVec3d& tau_segs) const;
};

}  // namespace sofa::component::cosserat::engine
