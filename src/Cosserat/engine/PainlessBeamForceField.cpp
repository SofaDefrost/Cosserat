#include <Cosserat/engine/PainlessBeamForceField.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/logging/Messaging.h>

#include <algorithm>   // std::max_element
#include <cmath>       // std::sqrt
#include <iomanip>     // std::setw, std::fixed, std::setprecision
#include <sstream>

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
                           "OUTPUT — differential torques on N segments (body frame)")),
      d_printEvery(initData(&d_printEvery, 100, "printEvery",
                            "Print debug info every N force computations. 0 = disabled.")),
      d_printPerSegment(initData(&d_printPerSegment, false, "printPerSegment",
                                 "Also print per-segment Gamma/Omega/f/tau when printing summary.")) {}

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

    // ── Stability estimate ────────────────────────────────────────────────────
    // Warn if the coupling stiffness (GA·h) will produce angular frequencies
    // above the Nyquist frequency for the current dt.  This is the known cause
    // of the zigzag instability in the Python explicit Euler integrator.
    {
        const auto& h_vec = l_state.get()->getRestLengths();
        if (!h_vec.empty()) {
            const double h0 = h_vec[0];
            const double GA = d_GA.getValue();
            const double GJ = d_GJ.getValue();
            const double EIy = d_EIy.getValue();
            // Effective angular frequency from coupling (GA shear → tau_lin path)
            // ω² ≈ GA·h / I_seg, but we don't have I_seg here.
            // Print the stiffness numbers and let the user judge.
            const double omega_ref_lin  = std::sqrt(GA / h0);   // [rad/s · sqrt(1/kg·m)]
            const double omega_ref_ang  = std::sqrt(EIy / (h0 * h0)); // [rad/s · sqrt(1/kg·m²)]
            msg_info() << "┌─ PainlessBeamForceField INIT ─────────────────────────────";
            msg_info() << "│  N=" << N << " segments   h=" << h0 << " m";
            msg_info() << "│  K_L : EA=" << d_EA.getValue() << " N   GA=" << GA << " N";
            msg_info() << "│  K_A : GJ=" << GJ << " N·m²   EIy=" << EIy
                       << " N·m²   EIz=" << d_EIz.getValue() << " N·m²";
            msg_info() << "│  sqrt(GA/h)  = " << omega_ref_lin
                       << "  [~angular freq if I_seg=1 kg·m²]";
            msg_info() << "│  sqrt(EIy/h²)= " << omega_ref_ang
                       << "  [~angular freq if I_seg=1 kg·m²]";
            msg_info() << "│  Stability check: for each I_seg, verify dt < 2/sqrt(GA*h/I_seg)";
            msg_info() << "│  printEvery=" << d_printEvery.getValue()
                       << "  printPerSegment=" << d_printPerSegment.getValue();
            msg_info() << "└───────────────────────────────────────────────────────────";
        }
    }

    // Enable automatic force update every animation step (for Python integrator)
    this->f_listening.setValue(true);
}

void PainlessBeamForceField::reinit() { init(); }

void PainlessBeamForceField::handleEvent(sofa::core::objectmodel::Event* e) {
    if (sofa::simulation::AnimateBeginEvent::checkEventType(e)) {
        // m_step is incremented inside computeForcesAndTorques.
        // (m_step+1) % freq == 0  → next increment will trigger the summary print.
        const int freq = d_printEvery.getValue();
        if (freq > 0 && ((m_step + 1) % freq == 0))
            msg_info() << "[handleEvent] AnimateBeginEvent — step=" << (m_step + 1)
                       << " → computeAndStoreForces() (summary will follow)";
        computeAndStoreForces();
    }
}

// ── Internal: stiffness matrices ──────────────────────────────────────────────

PainlessBeamForceField::Mat3x3d PainlessBeamForceField::buildK_L() const {
    Mat3x3d K;  K.clear();
    K[0][0] = d_EA.getValue();
    K[1][1] = d_GA.getValue();
    K[2][2] = d_GA.getValue();
    return K;
}

PainlessBeamForceField::Mat3x3d PainlessBeamForceField::buildK_A() const {
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

    // ── Per-step debug control ────────────────────────────────────────────────
    ++m_step;
    const int freq   = d_printEvery.getValue();
    const bool doPrint    = (freq > 0) && (m_step % freq == 0);
    const bool doPerSeg   = doPrint && d_printPerSegment.getValue();

    // Accumulators for summary (always computed, printed only if doPrint)
    double max_Gamma_norm = 0.0, max_Omega_norm = 0.0;
    double max_f_world    = 0.0, max_tau_lin    = 0.0;
    double max_tau_ang    = 0.0;
    double energy_lin = 0.0,    energy_ang = 0.0;

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

        // ── Torque on segment i from linear strain (coupling positions → SO3) ──
        //
        //   τ_linear_i = dx_body × f_body    (body frame of segment i)
        //
        // Derived from  δE_i = ξ · (dx_body × f_body)  under right perturbation
        // R_i → R_i exp(ε ξ).  This term drives R_i to follow the deformed beam
        // direction; without it SO3 DOFs are decoupled from position dynamics and
        // bending stiffness (EI/GJ) never activates.
        const SO3::Vector tau_lin = dx_body.cross(f_body);
        tau_segs[i] += Vec3d(tau_lin.x(), tau_lin.y(), tau_lin.z());

        const double e_lin_i = 0.5 * h[i] * (Gamma_i * KL_Gamma);
        energy_lin += e_lin_i;
        energy     += e_lin_i;

        // Accumulate extrema
        max_Gamma_norm = std::max(max_Gamma_norm, Gamma_i.norm());
        max_f_world    = std::max(max_f_world,    f_world.norm());
        max_tau_lin    = std::max(max_tau_lin,    Vec3d(tau_lin.x(), tau_lin.y(), tau_lin.z()).norm());

        if (doPerSeg) {
            // Rotation vector of R_i
            const auto omega_i = R[i].log();
            msg_info() << std::fixed << std::setprecision(4)
                       << "  [lin seg " << i << "]"
                       << "  h=" << h[i]
                       << "  dx=(" << dx.x() << "," << dx.y() << "," << dx.z() << ")"
                       << "  |dx|=" << dx.norm()
                       << "  dx_body=(" << dx_body.x() << "," << dx_body.y() << "," << dx_body.z() << ")"
                       << "  Gamma=(" << Gamma_i.x() << "," << Gamma_i.y() << "," << Gamma_i.z()
                       << ")  |Gamma|=" << Gamma_i.norm()
                       << "  f_body=(" << f_body.x() << "," << f_body.y() << "," << f_body.z() << ")"
                       << "  f_world=(" << f_world.x() << "," << f_world.y() << "," << f_world.z()
                       << ")  |f_world|=" << f_world.norm()
                       << "  tau_lin=(" << tau_lin.x() << "," << tau_lin.y() << "," << tau_lin.z()
                       << ")  |tau_lin|=" << Vec3d(tau_lin.x(), tau_lin.y(), tau_lin.z()).norm()
                       << "  omega_i=(" << omega_i.x() << "," << omega_i.y() << "," << omega_i.z() << ")"
                       << "  E_lin=" << e_lin_i;
        }
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
        const Vec3d tau_i = Jr_inv_neg * KA_Omega;
        tau_segs[i] -= tau_i;

        // τ on R_{i-1}
        const Mat3x3d Jr_inv_pos = CosseratIntrinsicState::getInverseLieJacobian(phi);
        const Vec3d tau_im1 = Jr_inv_pos * KA_Omega;
        tau_segs[i - 1] += tau_im1;

        const double e_ang_i = 0.5 * h_dual * (Omega * KA_Omega);
        energy_ang += e_ang_i;
        energy     += e_ang_i;

        max_Omega_norm = std::max(max_Omega_norm, Omega.norm());
        max_tau_ang    = std::max(max_tau_ang,    std::max(tau_i.norm(), tau_im1.norm()));

        if (doPerSeg) {
            msg_info() << std::fixed << std::setprecision(4)
                       << "  [ang node " << i << "]"
                       << "  h_dual=" << h_dual
                       << "  phi=(" << phi.x() << "," << phi.y() << "," << phi.z()
                       << ")  |phi|=" << phi.norm()
                       << "  Omega=(" << Omega.x() << "," << Omega.y() << "," << Omega.z()
                       << ")  |Omega|=" << Omega.norm()
                       << "  KA_Omega=(" << KA_Omega.x() << "," << KA_Omega.y() << "," << KA_Omega.z() << ")"
                       << "  tau_i=(" << tau_i.x() << "," << tau_i.y() << "," << tau_i.z()
                       << ")  |tau_i|=" << tau_i.norm()
                       << "  tau_im1=(" << tau_im1.x() << "," << tau_im1.y() << "," << tau_im1.z()
                       << ")  |tau_im1|=" << tau_im1.norm()
                       << "  E_ang=" << e_ang_i;
        }
    }

    // ── Summary print (every d_printEvery steps) ─────────────────────────────
    if (doPrint) {
        // Tip position
        const Vec3d& tip = pos[Np1 - 1];

        // Max accumulated torque per segment (lin + ang combined)
        double max_tau_total = 0.0;
        for (size_t i = 0; i < N; ++i)
            max_tau_total = std::max(max_tau_total, tau_segs[i].norm());

        // Max nodal force
        double max_f_node = 0.0;
        for (size_t j = 0; j < Np1; ++j)
            max_f_node = std::max(max_f_node, f_nodes[j].norm());

        // Rotation norms for all segments
        double max_omega_seg = 0.0;
        for (size_t i = 0; i < N; ++i) {
            const auto& omegas = state->d_orientations.getValue();
            max_omega_seg = std::max(max_omega_seg, omegas[i].norm());
        }

        msg_info() << "┌─ PainlessBeamForceField step=" << m_step << " ──────────────────────────────────";
        msg_info() << "│  tip pos : (" << std::fixed << std::setprecision(5)
                   << tip.x() << ", " << tip.y() << ", " << tip.z() << ")";
        msg_info() << "│  tip deflection |tip - tip0| = "
                   << (tip - pos[0]).norm() << " m";
        msg_info() << "│  Energy total=" << std::scientific << std::setprecision(3)
                   << energy << "  lin=" << energy_lin << "  ang=" << energy_ang;
        msg_info() << "│  max |Gamma|  (lin strain) = " << max_Gamma_norm
                   << "  [should be 0 at rest]";
        msg_info() << "│  max |Omega|  (ang strain) = " << max_Omega_norm
                   << "  [should be 0 at rest]";
        msg_info() << "│  max |f_world|(from lin)  = " << max_f_world << " N";
        msg_info() << "│  max |tau_lin|(coupling)  = " << max_tau_lin << " N·m"
                   << "  ← if >> max_tau_ang → instability risk";
        msg_info() << "│  max |tau_ang|(from ang)  = " << max_tau_ang << " N·m";
        msg_info() << "│  max |tau_total| (per seg) = " << max_tau_total << " N·m";
        msg_info() << "│  max |f_node|  (on nodes) = " << max_f_node << " N";
        msg_info() << "│  max |omega_seg| (rot vec) = " << max_omega_seg << " rad"
                   << "  [identity=0, warning if >π/2≈1.57]";
        msg_info() << "│  tau_lin / tau_ang ratio  = "
                   << (max_tau_ang > 1e-30 ? max_tau_lin / max_tau_ang : 0.0)
                   << "  [large → lin coupling dominates → instability]";
        msg_info() << "└──────────────────────────────────────────────────────────────";
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
        const auto rot_e = R[i].matrix();   // Eigen::Matrix3d
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
    sofa::core::MultiVecDerivId        /*fId*/) {

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

void PainlessBeamForceField::computeAndStoreForces() {
    if (!l_state.get()) {
        msg_warning() << "[computeAndStoreForces] l_state is null — skipping.";
        return;
    }
    VecVec3d& f_nodes  = *d_nodalForces.beginEdit();
    VecVec3d& tau_segs = *d_segmentTorques.beginEdit();
    computeForcesAndTorques(f_nodes, tau_segs);
    d_nodalForces.endEdit();
    d_segmentTorques.endEdit();
}

void PainlessBeamForceField::addDForce(
    const sofa::core::MechanicalParams* mparams,
    sofa::core::MultiVecDerivId        /*dfId*/) {

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
    const sofa::core::MechanicalParams* /*mparams*/,
    const sofa::core::behavior::MultiMatrixAccessor* /*matrix*/) {
    // Stub for painless/python-explicit: time integration is handled by a Python
    // explicit-Euler controller that calls computeDForces() directly.
    // Full sparse assembly is deferred to the painless/sofa-implicit branch.
}

SReal PainlessBeamForceField::getPotentialEnergy(
    const sofa::core::MechanicalParams* /*mparams*/) const {
    if (!l_state.get()) return 0.0;
    VecVec3d f_tmp, tau_tmp;
    return computeForcesAndTorques(f_tmp, tau_tmp);
}

}  // namespace sofa::component::cosserat::engine
