#include <Cosserat/engine/PainlessBeamForceField.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/logging/Messaging.h>

// Registration (called from initCosserat.cpp)
namespace Cosserat {
void registerPainlessBeamForceField(sofa::core::ObjectFactory* factory) {
    factory->registerObjects(
        sofa::core::ObjectRegistrationData("PainlessBeamForceField")
            .add<sofa::component::cosserat::engine::PainlessBeamForceField>());
}
}  // namespace Cosserat

namespace sofa::component::cosserat::engine {

// ── Constructor ───────────────────────────────────────────────────────────────

PainlessBeamForceField::PainlessBeamForceField()
    : d_EA(initData(&d_EA, 1.0e6, "EA",
                    "Axial stiffness E·A  [N]")),
      d_GA(initData(&d_GA, 1.0e5, "GA",
                    "Shear stiffness G·A  [N]  (same for both transverse directions)")),
      d_GJ(initData(&d_GJ, 1.0e5, "GJ",
                    "Torsion stiffness G·J  [N·m²]")),
      d_EIy(initData(&d_EIy, 1.0e4, "EIy",
                     "Bending stiffness E·Iy  [N·m²]  (y-axis)")),
      d_EIz(initData(&d_EIz, 1.0e4, "EIz",
                     "Bending stiffness E·Iz  [N·m²]  (z-axis)")),
      l_state(initLink("state",
                       "Link to the CosseratIntrinsicState that holds the beam DOFs")),
      d_nodalForces(initData(&d_nodalForces, "nodalForces",
                             "OUTPUT — elastic forces on N+1 position DOFs (world frame)")),
      d_segmentTorques(initData(&d_segmentTorques, "segmentTorques",
                                "OUTPUT — elastic torques on N SO3 DOFs (body frame)")) {}

// ── Lifecycle ─────────────────────────────────────────────────────────────────

void PainlessBeamForceField::init() {
    BaseForceField::init();

    if (!l_state.get()) {
        msg_error() << "No CosseratIntrinsicState linked via 'state'. "
                       "PainlessBeamForceField will produce no forces.";
        return;
    }

    const size_t N   = l_state.get()->getNbSegments();
    const size_t Np1 = l_state.get()->getPositions().size();

    if (Np1 != N + 1) {
        msg_error() << "CosseratIntrinsicState size mismatch: "
                    << "getNbSegments()=" << N << "  but positions.size()=" << Np1
                    << "  (expected N+1=" << N + 1 << ").";
        return;
    }

    msg_info() << "PainlessBeamForceField ready — N=" << N << " segments, "
               << "K_L=diag(" << d_EA.getValue() << ", " << d_GA.getValue() << ", "
               << d_GA.getValue() << ")  "
               << "K_A=diag(" << d_GJ.getValue() << ", " << d_EIy.getValue() << ", "
               << d_EIz.getValue() << ")";
}

void PainlessBeamForceField::reinit() { init(); }

// ── Internal: stiffness matrices ──────────────────────────────────────────────

Mat3x3d PainlessBeamForceField::buildK_L() const {
    Mat3x3d K;
    K.clear();
    K[0][0] = d_EA.getValue();
    K[1][1] = d_GA.getValue();
    K[2][2] = d_GA.getValue();
    return K;
}

Mat3x3d PainlessBeamForceField::buildK_A() const {
    Mat3x3d K;
    K.clear();
    K[0][0] = d_GJ.getValue();
    K[1][1] = d_EIy.getValue();
    K[2][2] = d_EIz.getValue();
    return K;
}

// ── Core computation ──────────────────────────────────────────────────────────

double PainlessBeamForceField::computeForcesAndTorques(
    VecVec3d& f_nodes,
    VecVec3d& tau_segs) const {
    const auto* state = l_state.get();
    if (!state) return 0.0;

    const auto& pos = state->getPositions();
    const auto& R   = state->getOrientations();
    const auto& h   = state->getRestLengths();

    const size_t N   = R.size();
    const size_t Np1 = pos.size();

    if (Np1 != N + 1 || h.size() != N) {
        msg_error() << "Size mismatch in computeForcesAndTorques — "
                    << "pos=" << Np1 << "  R=" << N << "  h=" << h.size();
        return 0.0;
    }

    f_nodes.assign(Np1, Vec3d(0.0, 0.0, 0.0));
    tau_segs.assign(N,  Vec3d(0.0, 0.0, 0.0));

    const Mat3x3d K_L = buildK_L();
    const Mat3x3d K_A = buildK_A();
    double energy = 0.0;

    // ── Linear strain (stretch + shear) ──────────────────────────────────────
    //
    //   Γ_i = R_i⁻¹ · (x_{i+1} − x_i) / h_i − e₁       (body frame)
    //   E_i = h_i/2 · Γ_i^T · K_L · Γ_i
    //
    //   Force on x_i     +=  R_i · K_L · Γ_i   (world frame)
    //   Force on x_{i+1} += −R_i · K_L · Γ_i
    //
    for (size_t i = 0; i < N; ++i) {
        if (h[i] < 1e-12) continue;

        // (x_{i+1} − x_i) rotated into segment body frame
        const Vec3d dx = pos[i + 1] - pos[i];
        const SO3::Vector dx_body =
            R[i].inverse().act(SO3::Vector(dx.x(), dx.y(), dx.z()));

        // Linear strain in body frame: Γ_i = dx_body/h_i − e₁
        const Vec3d Gamma_i(dx_body.x() / h[i] - 1.0,
                            dx_body.y() / h[i],
                            dx_body.z() / h[i]);

        // Elastic force in body frame, then rotated to world frame
        const Vec3d KL_Gamma = K_L * Gamma_i;
        const SO3::Vector f_body(KL_Gamma.x(), KL_Gamma.y(), KL_Gamma.z());
        const SO3::Vector f_world_eigen = R[i].act(f_body);
        const Vec3d f_world(f_world_eigen.x(), f_world_eigen.y(), f_world_eigen.z());

        f_nodes[i]     += f_world;   // left node: pulled toward right
        f_nodes[i + 1] -= f_world;   // right node: Newton's 3rd law

        energy += 0.5 * h[i] * (Gamma_i * KL_Gamma);   // dot product
    }

    // ── Angular strain (bending + torsion) ───────────────────────────────────
    //
    //   φ_i   = log(R_{i-1}^T · R_i)           (relative rotation, unscaled)
    //   h̃_i   = (h_{i-1} + h_i) / 2            (dual-edge length at node i)
    //   Ω_i   = φ_i / h̃_i                       (angular strain)
    //   E_i   = h̃_i/2 · Ω_i^T · K_A · Ω_i
    //
    //   Using J_r^{-T}(φ) = J_r^{-1}(−φ)  (identity in SO(3)):
    //
    //   τ on R_i     += −J_r^{-1}(−φ_i) · K_A · Ω_i
    //   τ on R_{i-1} += +J_r^{-1}( φ_i) · K_A · Ω_i
    //
    //   Torques are in the body frame of the respective segment.
    //   Boundary nodes i=0 and i=N are left at zero (BCs enforced externally).
    //
    for (size_t i = 1; i < N; ++i) {
        const double h_dual = (h[i - 1] + h[i]) * 0.5;
        if (h_dual < 1e-12) continue;

        // φ_i = log(R_{i-1}^T · R_i)
        const SO3 rel_R = R[i - 1].inverse() * R[i];
        const SO3::TangentVector phi_eigen = rel_R.log();
        const Vec3d phi(phi_eigen.x(), phi_eigen.y(), phi_eigen.z());

        // Angular strain Ω_i = φ_i / h̃_i
        const Vec3d Omega = phi / h_dual;

        // K_A · Ω_i  (diagonal K_A → element-wise product)
        const Vec3d KA_Omega = K_A * Omega;

        // Torque on R_i: −J_r^{-1}(−φ) · K_A · Ω
        const Mat3x3d Jr_inv_neg_phi =
            CosseratIntrinsicState::getInverseLieJacobian(-phi);
        tau_segs[i] -= Jr_inv_neg_phi * KA_Omega;

        // Torque on R_{i-1}: +J_r^{-1}(φ) · K_A · Ω
        const Mat3x3d Jr_inv_phi =
            CosseratIntrinsicState::getInverseLieJacobian(phi);
        tau_segs[i - 1] += Jr_inv_phi * KA_Omega;

        energy += 0.5 * h_dual * (Omega * KA_Omega);   // dot product
    }

    return energy;
}

// ── BaseForceField overrides ──────────────────────────────────────────────────

void PainlessBeamForceField::addForce(
    const sofa::core::MechanicalParams* /*mparams*/,
    sofa::core::MultiVecDerivId /*f*/,
    sofa::core::ConstMultiVecCoordId /*x*/,
    sofa::core::ConstMultiVecDerivId /*v*/) {
    if (!l_state.get()) return;

    VecVec3d& f_nodes  = *d_nodalForces.beginEdit();
    VecVec3d& tau_segs = *d_segmentTorques.beginEdit();

    computeForcesAndTorques(f_nodes, tau_segs);

    d_nodalForces.endEdit();
    d_segmentTorques.endEdit();

    // TODO(solver-integration): write f_nodes and tau_segs into the global
    // MultiVecDerivId 'f'.  Requires CosseratIntrinsicState to define the
    // standard write(VecDerivId) path used by SOFA implicit solvers.
    // Until then, forces are accessible via d_nodalForces / d_segmentTorques.
}

void PainlessBeamForceField::addDForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::core::MultiVecDerivId /*df*/,
    sofa::core::ConstMultiVecDerivId /*dx*/) {
    SOFA_UNUSED(mparams);

    // TODO(addDForce-linear): implement df = −kFactor · K_lin · dx for the
    // linear (stretch/shear) contribution.  The tangent block for segment i is:
    //
    //   K_lin_i = R_i · K_L · R_i^T / h_i   (3×3, world frame)
    //
    //   df(x_i)     += −kFactor · K_lin_i · (dx_{i+1} − dx_i)
    //   df(x_{i+1}) +=  kFactor · K_lin_i · (dx_{i+1} − dx_i)
    //
    // TODO(addDForce-angular): add the angular stiffness contribution once the
    // linear path has been validated.
    //
    // Blocked by: same VecDeriv access issue as addForce.
}

void PainlessBeamForceField::addKToMatrix(
    const sofa::core::MechanicalParams* mparams,
    const sofa::core::behavior::MultiMatrixAccessor* matrix) {
    if (!l_state.get()) return;

    // Attempt to retrieve the global matrix row/col offset for l_state.
    // This will only work once CosseratIntrinsicState is fully registered
    // as a MechanicalState with the SOFA matrix assembly pipeline.
    const sofa::core::behavior::MultiMatrixAccessor::MatrixRef mref =
        matrix->getMatrix(l_state.get());

    if (!mref.matrix) {
        // State is not yet wired into the matrix system.  Silently skip.
        return;
    }

    using BaseMatrix = sofa::linearalgebra::BaseMatrix;
    BaseMatrix* mat      = mref.matrix;
    const unsigned off   = mref.offset;
    const double kFactor =
        mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const auto& R = l_state.get()->getOrientations();
    const auto& h = l_state.get()->getRestLengths();
    const size_t N = R.size();

    if (h.size() != N) return;

    const Mat3x3d K_L = buildK_L();

    // ── Linear stiffness blocks ───────────────────────────────────────────────
    //
    //   For segment i (i = 0..N-1):
    //
    //     K_world_i = R_i · K_L · R_i^T / h_i    (3×3 block, world frame)
    //
    //   Contributions to the global 3(N+1) × 3(N+1) position sub-matrix:
    //
    //     K[3i  : 3i+3,  3i  : 3i+3]   += +kFactor · K_world_i
    //     K[3(i+1):…,  3(i+1):…]       += +kFactor · K_world_i
    //     K[3i  : 3i+3, 3(i+1):…]      += -kFactor · K_world_i
    //     K[3(i+1):…,  3i  : 3i+3]     += -kFactor · K_world_i   (symmetry)
    //
    for (size_t i = 0; i < N; ++i) {
        if (h[i] < 1e-12) continue;

        // K_world_i = R_i · K_L · R_i^T / h_i
        // Build R_i as a 3×3 SOFA matrix
        const auto rot_mat = R[i].toRotationMatrix();  // Eigen::Matrix3d
        Mat3x3d R_sofa;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                R_sofa[r][c] = rot_mat(r, c);

        // R · K_L
        const Mat3x3d R_KL = R_sofa * K_L;
        // (R · K_L) · R^T = R · K_L · R^T
        Mat3x3d K_world;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                double sum = 0.0;
                for (int k = 0; k < 3; ++k)
                    sum += R_KL[r][k] * R_sofa[c][k];  // R_sofa[c][k] = R^T[k][c]
                K_world[r][c] = sum / h[i];
            }

        const unsigned row_i   = off + static_cast<unsigned>(3 * i);
        const unsigned row_ip1 = off + static_cast<unsigned>(3 * (i + 1));

        for (unsigned r = 0; r < 3; ++r)
            for (unsigned c = 0; c < 3; ++c) {
                const double val = kFactor * K_world[r][c];
                // diagonal blocks
                mat->add(row_i   + r, row_i   + c, +val);
                mat->add(row_ip1 + r, row_ip1 + c, +val);
                // off-diagonal blocks (negative coupling)
                mat->add(row_i   + r, row_ip1 + c, -val);
                mat->add(row_ip1 + r, row_i   + c, -val);
            }
    }

    // TODO(addKToMatrix-angular): add the 3×3 angular stiffness blocks
    //   K_ang_i = J_r^{-T}(φ_i) · K_A · J_r^{-1}(φ_i) / h̃_i
    // into the SO3 DOF sub-matrix once the DOF index mapping for orientations
    // is established in the matrix assembly pipeline.
}

SReal PainlessBeamForceField::getPotentialEnergy(
    const sofa::core::MechanicalParams* /*mparams*/,
    sofa::core::ConstMultiVecCoordId /*x*/) const {
    if (!l_state.get()) return 0.0;

    // computeForcesAndTorques is logically const here (output buffers are
    // temporary).  We use local temporaries rather than modifying Data.
    VecVec3d f_tmp, tau_tmp;
    return computeForcesAndTorques(f_tmp, tau_tmp);
}

}  // namespace sofa::component::cosserat::engine
