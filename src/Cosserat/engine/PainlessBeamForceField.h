#pragma once

#include <Cosserat/engine/CosseratIntrinsicState.h>
#include <liegroups/SO3.h>

#include <sofa/core/behavior/BaseForceField.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Link.h>
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
 * ## Tangent stiffness (addDForce)
 *
 * ### Linear material stiffness (per segment i):
 * @code
 *   K_world_i = R_i · K_L · R_i^T / h_i
 *   df(x_i)     += +kF · K_world_i · (dx_{i+1} − dx_i)
 *   df(x_{i+1}) += −kF · K_world_i · (dx_{i+1} − dx_i)
 * @endcode
 *
 * ### Angular material stiffness (per interior node i), with
 *   A = J_r^{-1}(−φ_i), B = J_r^{-1}(φ_i):
 * @code
 *   dτ(R_i)     += kF · [−(A·K_A·B/h̃_i)·dω_i + (A·K_A·A/h̃_i)·dω_{i−1}]
 *   dτ(R_{i-1}) += kF · [+(B·K_A·B/h̃_i)·dω_i − (B·K_A·A/h̃_i)·dω_{i−1}]
 * @endcode
 * (material stiffness only; geometric stiffness terms from ∂J_r^{-1}/∂φ are TODO)
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
    sofa::core::objectmodel::Data<double> d_GA;   ///< Shear stiffness  G·A
    sofa::core::objectmodel::Data<double> d_GJ;   ///< Torsion stiffness G·J
    sofa::core::objectmodel::Data<double> d_EIy;  ///< Bending stiffness E·Iy
    sofa::core::objectmodel::Data<double> d_EIz;  ///< Bending stiffness E·Iz

    // ── Link to the mechanical state ───────────────────────────────────────────
    sofa::core::objectmodel::SingleLink<
        PainlessBeamForceField,
        CosseratIntrinsicState,
        sofa::core::objectmodel::BaseLink::FLAG_STOREPATH |
        sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK
    > l_state;

    // ── Force/torque outputs ───────────────────────────────────────────────────
    /// Elastic forces on N+1 position DOFs (world frame).
    sofa::core::objectmodel::Data<VecVec3d> d_nodalForces;
    /// Elastic torques on N SO3 DOFs (body frame).
    sofa::core::objectmodel::Data<VecVec3d> d_segmentTorques;

    // ── Differential force inputs/outputs (for Python integration and addDForce)
    /// INPUT  — position displacements Δx_i ∈ R³ for i = 0..N  (world frame).
    sofa::core::objectmodel::Data<VecVec3d> d_dx_positions;
    /// INPUT  — angular displacements Δω_i ∈ so(3) for i = 0..N-1 (body frame).
    sofa::core::objectmodel::Data<VecVec3d> d_dx_angles;
    /// OUTPUT — differential forces  df on N+1 position DOFs (world frame).
    sofa::core::objectmodel::Data<VecVec3d> d_df_positions;
    /// OUTPUT — differential torques dτ on N SO3 DOFs (body frame).
    sofa::core::objectmodel::Data<VecVec3d> d_df_angles;

    // ─────────────────────────────────────────────────────────────────────────
    PainlessBeamForceField();
    ~PainlessBeamForceField() override = default;

    void init() override;
    void reinit() override;

    // ── BaseForceField interface ───────────────────────────────────────────────
    //
    // Signatures match sofa::core::behavior::BaseForceField (SOFA v24+):
    //   addForce(mparams, fId)        — state read/write via MultiVecDerivId
    //   addDForce(mparams, dfId)      — differential force
    //   getPotentialEnergy(mparams)   — scalar energy
    //   addKToMatrix(mparams, matrix) — sparse stiffness assembly
    //
    void addForce(const sofa::core::MechanicalParams* mparams,
                  sofa::core::MultiVecDerivId fId) override;

    void addDForce(const sofa::core::MechanicalParams* mparams,
                   sofa::core::MultiVecDerivId dfId) override;

    void addKToMatrix(const sofa::core::MechanicalParams* mparams,
                      const sofa::core::behavior::MultiMatrixAccessor* matrix) override;

    SReal getPotentialEnergy(const sofa::core::MechanicalParams* mparams) const override;

    // ── Direct Python/C++ API ──────────────────────────────────────────────────

    /**
     * @brief Compute differential forces from explicit displacement vectors.
     *
     * Reads d_dx_positions and d_dx_angles, writes d_df_positions and
     * d_df_angles.  Exposed to Python via the Cosserat pybind11 module so that
     * a Python explicit-Euler integrator can call it directly.
     *
     * @param kFactor  Stiffness scale factor (1.0 for explicit Euler).
     */
    void computeDForcesFromData(double kFactor = 1.0);

    /**
     * @brief Core differential force computation (standalone, no Data I/O).
     *
     * @param dx_pos   Position displacements  Δx_i  (size N+1).
     * @param dx_ang   Angular displacements   Δω_i  (size N, so3 body frame).
     * @param kFactor  Stiffness scale factor.
     * @param[out] df_pos  Differential forces  on N+1 nodes (world frame).
     * @param[out] df_ang  Differential torques on N segments (body frame).
     */
    void computeDForces(const VecVec3d& dx_pos,
                        const VecVec3d& dx_ang,
                        double          kFactor,
                        VecVec3d&       df_pos,
                        VecVec3d&       df_ang) const;

   private:
    Mat3x3d buildK_L() const;
    Mat3x3d buildK_A() const;

    double computeForcesAndTorques(VecVec3d& f_nodes,
                                   VecVec3d& tau_segs) const;
};

}  // namespace sofa::component::cosserat::engine
