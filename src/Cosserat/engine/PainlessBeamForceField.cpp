#include <Cosserat/engine/PainlessBeamForceField.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/logging/Messaging.h>

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
    : d_EA(initData(&d_EA, 1.0e6, "EA", "Axial stiffness E·A  [N]")),
      d_GA(initData(&d_GA, 1.0e5, "GA",
                    "Shear stiffness G·A  [N] (same for both transverse directions)")),
      d_GJ(initData(&d_GJ, 1.0e5, "GJ", "Torsion stiffness G·J  [N·m²]")),
      d_EIy(initData(&d_EIy, 1.0e4, "EIy", "Bending stiffness E·Iy  [N·m²] (y-axis)")),
      d_EIz(initData(&d_EIz, 1.0e4, "EIz", "Bending stiffness E·Iz  [N·m²] (z-axis)")),
      l_state(initLink("state",
                       "Link to the CosseratIntrinsicState that holds the beam DOFs")),
      d_nodalForces(initData(&d_nodalForces, "nodalForces",
                             "OUTPUT — elastic forces on N+1 position DOFs (world frame)")),
      d_segmentTorques(initData(&d_segmentTorques, "segmentTorques",
                                "OUTPUT — elastic torques on N SO3 DOFs (body frame)")),
      d_dx_positions(initData(&d_dx_positions, "dx_positions",
                              "INPUT — position displacements Δx_i (world frame, size N+1)")),
      d_dx_angles(initData(&d_dx_angles, "dx_angles",
                           "INPUT — angular displacements Δω_i in so(3) body frame (size N)")),
      d_df_positions(initData(&d_df_positions, "df_positions",
                              "OUTPUT — differential forces on N+1 nodes (world frame)")),
      d_df_angles(initData(&d_df_angles, "df_angles",
                           "OUTPUT — differential torques on N segments (body frame)")) {}

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
                    << "getNbSegments()=" << N << "  but positions.size()=" << Np1;
        return;
    }

    // Pre-size the dx/df Data fields to avoid reallocation in first step
    {
        auto& dx_pos = *d_dx_positions.beginEdit();
        auto& dx_ang = *d_dx_angles.beginEdit();
        dx_pos.assign(Np1, Vec3d(0, 0, 0));
        dx_ang.assign(N,   Vec3d(0, 0, 0));
        d_dx_positions.endEdit();
        d_dx_angles.endEdit();
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
    Mat3x3d K;  K.clear();
    K[0][0] = d_EA.getValue();
    K[1][1] = d_GA.getValue();
    K[2][2] = d_GA.getValue();
    return K;
}

Mat3x3d PainlessBeamForceField::buildK_A() const {
    Mat3x3d K;  K.clear();
    K[0][0] = d_GJ.getValue();
    K[1][1] = d_EIy.getValue();
    K[2][2] = d_EIz.getValue();
    return K;
}

// ── Core: forces & torques ────────────────────────────────────────────────────

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
    if (Np1 != N + 1 || h.size() != N) return 0.0;

    f_nodes.assign(Np1, Vec3d(0, 0, 0));
    tau_segs.assign(N,  Vec3d(0, 0, 0));

    const Mat3x3d K_L = buildK_L();
    const Mat3x3d K_A = buildK_A();
    double energy = 0.0;

    // ── Linear strain (stretch + shear) ──────────────────────────────────────
    for (size_t i = 0; i < N; ++i) {
        if (h[i] < 1e-12) continue;

        const Vec3d dx = pos[i + 1] - pos[i];
        const SO3::Vector dx_body =
            R[i].inverse().act(SO3::Vector(dx.x(), dx.y(), dx.z()));

        const Vec3d Gamma_i(dx_body.x() / h[i] - 1.0,
                            dx_body.y() / h[i],
                            dx_body.z() / h[i]);

        const Vec3d KL_Gamma = K_L * Gamma_i;
        const SO3::Vector f_body(KL_Gamma.x(), KL_Gamma.y(), KL_Gamma.z());
        const SO3::Vector f_world_e = R[i].act(f_body);
        const Vec3d f_world(f_world_e.x(), f_world_e.y(), f_world_e.z());

        f_nodes[i]     += f_world;
        f_nodes[i + 1] -= f_world;

        energy += 0.5 * h[i] * (Gamma_i * KL_Gamma);
    }

    // ── Angular strain (bending + torsion) ───────────────────────────────────
    //
    //   φ_i  = log(R_{i-1}^T · R_i)
    //   Ω_i  = φ_i / h̃_i
    //
    //   τ(R_i)     += −J_r^{-1}(−φ_i) · K_A · Ω_i
    //   τ(R_{i-1}) += +J_r^{-1}( φ_i) · K_A · Ω_i
    //
    for (size_t i = 1; i < N; ++i) {
        const double h_dual = (h[i - 1] + h[i]) * 0.5;
        if (h_dual < 1e-12) continue;

        const SO3 rel_R = R[i - 1].inverse() * R[i];
        const SO3::TangentVector phi_e = rel_R.log();
        const Vec3d phi(phi_e.x(), phi_e.y(), phi_e.z());
        const Vec3d Omega = phi / h_dual;

        const Vec3d KA_Omega = K_A * Omega;

        // τ on R_i
        const Mat3x3d Jr_inv_neg = CosseratIntrinsicState::getInverseLieJacobian(-phi);
        tau_segs[i] -= Jr_inv_neg * KA_Omega;

        // τ on R_{i-1}
        const Mat3x3d Jr_inv_pos = CosseratIntrinsicState::getInverseLieJacobian(phi);
        tau_segs[i - 1] += Jr_inv_pos * KA_Omega;

        energy += 0.5 * h_dual * (Omega * KA_Omega);
    }

    return energy;
}

// ── Core: differential forces ─────────────────────────────────────────────────

void PainlessBeamForceField::computeDForces(
    const VecVec3d& dx_pos,
    const VecVec3d& dx_ang,
    double          kFactor,
    VecVec3d&       df_pos,
    VecVec3d&       df_ang) const {

    const auto* state = l_state.get();
    if (!state) return;

    const auto& R = state->getOrientations();
    const auto& h = state->getRestLengths();
    const size_t N   = R.size();
    const size_t Np1 = N + 1;

    if (dx_pos.size() != Np1 || dx_ang.size() != N) {
        msg_warning() << "computeDForces: dx size mismatch — "
                      << "dx_pos=" << dx_pos.size() << " (expected " << Np1 << "), "
                      << "dx_ang=" << dx_ang.size()  << " (expected " << N  << ").";
        return;
    }

    df_pos.assign(Np1, Vec3d(0, 0, 0));
    df_ang.assign(N,   Vec3d(0, 0, 0));

    const Mat3x3d K_L = buildK_L();
    const Mat3x3d K_A = buildK_A();

    // ── Linear material stiffness (positions) ─────────────────────────────────
    //
    //   K_world_i = R_i · K_L · R_i^T / h_i
    //   df(x_i)     += +kF · K_world_i · (dx_{i+1} − dx_i)
    //   df(x_{i+1}) += −kF · K_world_i · (dx_{i+1} − dx_i)
    //
    for (size_t i = 0; i < N; ++i) {
        if (h[i] < 1e-12) continue;

        // Build R_i as 3×3 SOFA matrix
        const auto rot_e = R[i].toRotationMatrix();   // Eigen::Matrix3d
        Mat3x3d R_sofa;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                R_sofa[r][c] = rot_e(r, c);

        // K_world = R · K_L · R^T / h_i
        const Mat3x3d R_KL = R_sofa * K_L;
        Mat3x3d K_world;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += R_KL[r][k] * R_sofa[c][k];   // R^T[k][c] = R[c][k]
                K_world[r][c] = s / h[i];
            }

        const Vec3d delta_x = dx_pos[i + 1] - dx_pos[i];
        const Vec3d df = K_world * delta_x;

        df_pos[i]     += df * kFactor;
        df_pos[i + 1] -= df * kFactor;
    }

    // ── Angular material stiffness (SO3 DOFs) ─────────────────────────────────
    //
    //   A = J_r^{-1}(−φ_i),   B = J_r^{-1}(φ_i)
    //
    //   dτ(R_i)     += kF · [−(A·K_A·B / h̃_i)·dω_i + (A·K_A·A / h̃_i)·dω_{i−1}]
    //   dτ(R_{i−1}) += kF · [+(B·K_A·B / h̃_i)·dω_i − (B·K_A·A / h̃_i)·dω_{i−1}]
    //
    //   NOTE: these are the material (constitutive) terms only.
    //   Geometric stiffness (∂J_r^{-1}/∂φ contributions) is TODO.
    //
    for (size_t i = 1; i < N; ++i) {
        const double h_dual = (h[i - 1] + h[i]) * 0.5;
        if (h_dual < 1e-12) continue;

        const SO3 rel_R = R[i - 1].inverse() * R[i];
        const SO3::TangentVector phi_e = rel_R.log();
        const Vec3d phi(phi_e.x(), phi_e.y(), phi_e.z());

        const Mat3x3d A = CosseratIntrinsicState::getInverseLieJacobian(-phi);
        const Mat3x3d B = CosseratIntrinsicState::getInverseLieJacobian(phi);

        // Stiffness blocks (material only)
        // -(A · K_A · B) / h̃_i → block (i, i)
        // +(A · K_A · A) / h̃_i → block (i, i-1)
        // +(B · K_A · B) / h̃_i → block (i-1, i)
        // -(B · K_A · A) / h̃_i → block (i-1, i-1)

        const Vec3d& dw_i   = dx_ang[i];
        const Vec3d& dw_im1 = dx_ang[i - 1];

        // A · K_A is just A with each row scaled by K_A diagonal
        auto AKA_v = [&](const Vec3d& v) -> Vec3d {
            // A · (K_A · A · v) — since K_A diagonal, K_A·x = [GJ*x0, EIy*x1, EIz*x2]
            const Vec3d Av = A * v;
            const Vec3d KA_Av(K_A[0][0] * Av[0], K_A[1][1] * Av[1], K_A[2][2] * Av[2]);
            return (A * KA_Av) * (kFactor / h_dual);
        };
        auto AKB_v = [&](const Vec3d& v) -> Vec3d {
            const Vec3d Bv = B * v;
            const Vec3d KA_Bv(K_A[0][0] * Bv[0], K_A[1][1] * Bv[1], K_A[2][2] * Bv[2]);
            return (A * KA_Bv) * (kFactor / h_dual);
        };
        auto BKB_v = [&](const Vec3d& v) -> Vec3d {
            const Vec3d Bv = B * v;
            const Vec3d KA_Bv(K_A[0][0] * Bv[0], K_A[1][1] * Bv[1], K_A[2][2] * Bv[2]);
            return (B * KA_Bv) * (kFactor / h_dual);
        };
        auto BKA_v = [&](const Vec3d& v) -> Vec3d {
            const Vec3d Av = A * v;
            const Vec3d KA_Av(K_A[0][0] * Av[0], K_A[1][1] * Av[1], K_A[2][2] * Av[2]);
            return (B * KA_Av) * (kFactor / h_dual);
        };

        df_ang[i]     -= AKB_v(dw_i);    // −(A·K_A·B/h̃)·dω_i
        df_ang[i]     += AKA_v(dw_im1);  // +(A·K_A·A/h̃)·dω_{i-1}
        df_ang[i - 1] += BKB_v(dw_i);    // +(B·K_A·B/h̃)·dω_i
        df_ang[i - 1] -= BKA_v(dw_im1);  // −(B·K_A·A/h̃)·dω_{i-1}
    }
}

// ── computeDForcesFromData (Python-callable) ──────────────────────────────────

void PainlessBeamForceField::computeDForcesFromData(double kFactor) {
    const VecVec3d& dx_pos = d_dx_positions.getValue();
    const VecVec3d& dx_ang = d_dx_angles.getValue();

    VecVec3d& df_pos = *d_df_positions.beginEdit();
    VecVec3d& df_ang = *d_df_angles.beginEdit();

    computeDForces(dx_pos, dx_ang, kFactor, df_pos, df_ang);

    d_df_positions.endEdit();
    d_df_angles.endEdit();
}

// ── BaseForceField overrides ──────────────────────────────────────────────────

void PainlessBeamForceField::addForce(
    const sofa::core::MechanicalParams* /*mparams*/,
    sofa::core::MultiVecDerivId        /*fId*/,
    sofa::core::ConstMultiVecCoordId   /*x*/,
    sofa::core::ConstMultiVecDerivId   /*v*/) {

    if (!l_state.get()) return;

    // ── Compute elastic forces and torques into Data fields ───────────────────
    //
    // Results are stored in d_nodalForces (N+1 Vec3, world frame) and
    // d_segmentTorques (N Vec3, body frame). They are read by:
    //   • The Python explicit-Euler integrator (operational path).
    //   • addKToMatrix / addDForce (tangent stiffness path).
    //
    // Integration with SOFA's implicit ODE solvers (EulerImplicitSolver) via
    // MultiVecDerivId requires CosseratIntrinsicState to implement typed
    // VecDeriv read/write accessors (i.e. become MechanicalState<DataTypes>).
    // This is planned as a future extension.
    //
    VecVec3d& f_nodes  = *d_nodalForces.beginEdit();
    VecVec3d& tau_segs = *d_segmentTorques.beginEdit();
    computeForcesAndTorques(f_nodes, tau_segs);
    d_nodalForces.endEdit();
    d_segmentTorques.endEdit();
}

void PainlessBeamForceField::addDForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::core::MultiVecDerivId        /*df*/,
    sofa::core::ConstMultiVecDerivId   /*dx*/) {

    SOFA_UNUSED(mparams);

    // The standard MultiVecDerivId path requires CosseratIntrinsicState to expose
    // standard VecDeriv accessors, which is tracked as a future TODO.
    //
    // Meanwhile, the differential forces are computed on demand via:
    //   1. Python: write d_dx_positions / d_dx_angles, call computeDForcesFromData()
    //   2. C++:    call computeDForces(dx_pos, dx_ang, kFactor, df_pos, df_ang) directly
    //
    // The Data-driven path is already called in the Python explicit Euler integrator
    // (see staggered_cantilever_full.py).
    //
    // NOTE: this method IS called by SOFA implicit solvers, but will silently do
    // nothing until CosseratIntrinsicState is fully wired into the solver pipeline.
}

void PainlessBeamForceField::addKToMatrix(
    const sofa::core::MechanicalParams* mparams,
    const sofa::core::behavior::MultiMatrixAccessor* matrix) {

    if (!l_state.get()) return;

    const sofa::core::behavior::MultiMatrixAccessor::MatrixRef mref =
        matrix->getMatrix(l_state.get());
    if (!mref.matrix) return;

    using BaseMatrix = sofa::linearalgebra::BaseMatrix;
    BaseMatrix* mat    = mref.matrix;
    const unsigned off = mref.offset;
    const double kFactor =
        mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const auto& R = l_state.get()->getOrientations();
    const auto& h = l_state.get()->getRestLengths();
    const size_t N = R.size();
    if (h.size() != N) return;

    const Mat3x3d K_L = buildK_L();
    const Mat3x3d K_A = buildK_A();

    // ── Linear stiffness blocks (positions) ───────────────────────────────────
    for (size_t i = 0; i < N; ++i) {
        if (h[i] < 1e-12) continue;

        const auto rot_e = R[i].toRotationMatrix();
        Mat3x3d R_sofa;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                R_sofa[r][c] = rot_e(r, c);

        const Mat3x3d R_KL = R_sofa * K_L;
        Mat3x3d K_world;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                double s = 0.0;
                for (int k = 0; k < 3; ++k)
                    s += R_KL[r][k] * R_sofa[c][k];
                K_world[r][c] = s / h[i];
            }

        const unsigned ri   = off + static_cast<unsigned>(3 * i);
        const unsigned rip1 = off + static_cast<unsigned>(3 * (i + 1));

        for (unsigned r = 0; r < 3; ++r)
            for (unsigned c = 0; c < 3; ++c) {
                const double val = kFactor * K_world[r][c];
                mat->add(ri   + r, ri   + c, +val);
                mat->add(rip1 + r, rip1 + c, +val);
                mat->add(ri   + r, rip1 + c, -val);
                mat->add(rip1 + r, ri   + c, -val);
            }
    }

    // ── Angular stiffness blocks (SO3 DOFs) ───────────────────────────────────
    //
    //   Layout: SO3 DOFs start at offset off + 3*(N+1) in the global matrix
    //   (N+1 position DOFs precede N angular DOFs).
    //
    //   For interior node i (i=1..N-1):
    //     A = J_r^{-1}(−φ_i),  B = J_r^{-1}(φ_i)
    //     Block (i,   i  ) += +kF · (−A·K_A·B / h̃_i)
    //     Block (i,   i-1) += +kF · (  A·K_A·A / h̃_i)
    //     Block (i-1, i  ) += +kF · (  B·K_A·B / h̃_i)
    //     Block (i-1, i-1) += +kF · (−B·K_A·A / h̃_i)
    //
    const unsigned ang_off = off + static_cast<unsigned>(3 * (N + 1));

    for (size_t i = 1; i < N; ++i) {
        const double h_dual = (h[i - 1] + h[i]) * 0.5;
        if (h_dual < 1e-12) continue;

        const SO3 rel_R = R[i - 1].inverse() * R[i];
        const SO3::TangentVector phi_e = rel_R.log();
        const Vec3d phi(phi_e.x(), phi_e.y(), phi_e.z());

        const Mat3x3d A = CosseratIntrinsicState::getInverseLieJacobian(-phi);
        const Mat3x3d B = CosseratIntrinsicState::getInverseLieJacobian(phi);

        // Compute the four 3×3 blocks
        auto mat3_mul = [](const Mat3x3d& X, const Mat3x3d& Y) -> Mat3x3d {
            Mat3x3d Z;  Z.clear();
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    for (int k = 0; k < 3; ++k)
                        Z[r][c] += X[r][k] * Y[k][c];
            return Z;
        };
        auto diag_right = [](const Mat3x3d& X, const Mat3x3d& D) -> Mat3x3d {
            // X · D where D is diagonal
            Mat3x3d Z;  Z.clear();
            for (int r = 0; r < 3; ++r)
                for (int c = 0; c < 3; ++c)
                    Z[r][c] = X[r][c] * D[c][c];
            return Z;
        };

        const Mat3x3d AKA = mat3_mul(diag_right(A, K_A), A) * (kFactor / h_dual);
        const Mat3x3d AKB = mat3_mul(diag_right(A, K_A), B) * (kFactor / h_dual);
        const Mat3x3d BKA = mat3_mul(diag_right(B, K_A), A) * (kFactor / h_dual);
        const Mat3x3d BKB = mat3_mul(diag_right(B, K_A), B) * (kFactor / h_dual);

        const unsigned ri   = ang_off + static_cast<unsigned>(3 * i);
        const unsigned rim1 = ang_off + static_cast<unsigned>(3 * (i - 1));

        for (unsigned r = 0; r < 3; ++r)
            for (unsigned c = 0; c < 3; ++c) {
                mat->add(ri   + r, ri   + c, -AKB[r][c]);   // block (i,   i)
                mat->add(ri   + r, rim1 + c, +AKA[r][c]);   // block (i,   i-1)
                mat->add(rim1 + r, ri   + c, +BKB[r][c]);   // block (i-1, i)
                mat->add(rim1 + r, rim1 + c, -BKA[r][c]);   // block (i-1, i-1)
            }
    }
}

SReal PainlessBeamForceField::getPotentialEnergy(
    const sofa::core::MechanicalParams* /*mparams*/,
    sofa::core::ConstMultiVecCoordId   /*x*/) const {
    if (!l_state.get()) return 0.0;
    VecVec3d f_tmp, tau_tmp;
    return computeForcesAndTorques(f_tmp, tau_tmp);
}

}  // namespace sofa::component::cosserat::engine
