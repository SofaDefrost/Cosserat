#include <Cosserat/engine/CosseratIntrinsicState.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/logging/Messaging.h>

#include <cmath>
#include <iomanip>
#include <sstream>

namespace Cosserat {
void registerCosseratIntrinsicState(sofa::core::ObjectFactory *factory) {
    factory->registerObjects(sofa::core::ObjectRegistrationData("CosseratIntrinsicState")
                                 .add<sofa::component::cosserat::engine::CosseratIntrinsicState>());
}
}  // namespace Cosserat

namespace sofa::component::cosserat::engine {

CosseratIntrinsicState::CosseratIntrinsicState()
    : d_positions(initData(&d_positions, "positions", "Positions of the N+1 nodes")),
      d_orientations(
          initData(&d_orientations, "orientations", "Orientations (SO3) of the N segments")),
      d_velocities(initData(&d_velocities, "velocities",
                            "Translational velocities of nodes in tangent space")),
      d_angularVelocities(initData(&d_angularVelocities, "angularVelocities",
                                   "Angular velocities of segments in tangent space")),
      d_accelerations(initData(&d_accelerations, "accelerations",
                               "Translational accelerations of nodes in tangent space")),
      d_angularAccelerations(initData(&d_angularAccelerations, "angularAccelerations",
                                      "Angular accelerations of segments in tangent space")) {}

sofa::type::vector<CosseratIntrinsicState::SO3>
CosseratIntrinsicState::getOrientations() const {
    const auto& omegas = d_orientations.getValue();
    sofa::type::vector<SO3> Rs(omegas.size());
    for (size_t i = 0; i < omegas.size(); ++i)
        Rs[i] = SO3::exp(SO3::TangentVector(omegas[i].x(), omegas[i].y(), omegas[i].z()));
    return Rs;
}

void CosseratIntrinsicState::init() {
    BaseObject::init();

    const auto& pos    = d_positions.getValue();
    const auto& omegas = d_orientations.getValue();
    size_t nbSections  = omegas.size();

    rest_lengths.clear();
    rest_lengths.reserve(nbSections);

    msg_info() << "┌─ CosseratIntrinsicState INIT ──────────────────────────────";
    msg_info() << "│  d_positions.size()    = " << pos.size()    << "  (expected N+1=" << nbSections+1 << ")";
    msg_info() << "│  d_orientations.size() = " << omegas.size() << "  (expected N="   << nbSections   << ")";

    bool length_ok = true;
    for (size_t i = 0; i < nbSections; ++i) {
        if (i + 1 < pos.size()) {
            double L = (pos[i+1] - pos[i]).norm();
            rest_lengths.push_back(L);
            if (L < 1e-12) {
                msg_warning() << "│  ⚠  rest_lengths[" << i << "] = " << L
                              << " (nearly zero — check positions)";
                length_ok = false;
            }
        } else {
            msg_warning() << "│  ⚠  Node index " << i+1 << " out of bounds — defaulting rest_length=1.0";
            rest_lengths.push_back(1.0);
            length_ok = false;
        }
    }

    // Print rest lengths summary
    double h_min = 1e30, h_max = 0.0, h_sum = 0.0;
    for (double h : rest_lengths) { h_min = std::min(h_min, h); h_max = std::max(h_max, h); h_sum += h; }
    msg_info() << "│  rest_lengths : min=" << h_min << "  max=" << h_max
               << "  total=" << h_sum << " m"
               << (length_ok ? "" : "  ← WARNINGS above");

    // Print initial node positions
    msg_info() << "│  Initial node positions (world frame):";
    for (size_t i = 0; i < pos.size(); ++i)
        msg_info() << "│    x[" << i << "] = ("
                   << std::fixed << std::setprecision(5)
                   << pos[i].x() << ", " << pos[i].y() << ", " << pos[i].z() << ")";

    // Print initial orientations (rotation vectors)
    msg_info() << "│  Initial orientations ω = log(R)  [identity = (0,0,0)]:";
    double max_omega0 = 0.0;
    for (size_t i = 0; i < omegas.size(); ++i) {
        msg_info() << "│    ω[" << i << "] = ("
                   << omegas[i].x() << ", " << omegas[i].y() << ", " << omegas[i].z()
                   << ")  |ω|=" << omegas[i].norm();
        max_omega0 = std::max(max_omega0, omegas[i].norm());
    }
    if (max_omega0 < 1e-10)
        msg_info() << "│  ✓ All orientations = identity (rotation vectors ≈ 0)";
    else
        msg_info() << "│  ⚠  Non-zero initial orientations (max |ω|=" << max_omega0 << ")";

    msg_info() << "└──────────────────────────────────────────────────────────────";
}

void CosseratIntrinsicState::updateStrainsCache() {
    bool dirty = false;
    const bool pos_changed = (d_positions.getCounter() != last_positions_counter);
    const bool ori_changed = (d_orientations.getCounter() != last_orientations_counter);

    if (pos_changed) {
        dirty = true;
        last_positions_counter = d_positions.getCounter();
    }
    if (ori_changed) {
        dirty = true;
        last_orientations_counter = d_orientations.getCounter();
    }

    if (!dirty) return;

    ++cache_update_count_;
    // Print the first 3 cache updates (init check) then every 500 to avoid spam
    const bool printCache = (cache_update_count_ <= 3) || (cache_update_count_ % 500 == 0);
    if (printCache) {
        msg_info() << "[updateStrainsCache] #" << cache_update_count_
                   << "  pos_changed=" << pos_changed
                   << "  ori_changed=" << ori_changed
                   << "  pos_counter=" << last_positions_counter
                   << "  ori_counter=" << last_orientations_counter;
    }

    const auto& pos    = d_positions.getValue();
    const auto& omegas = d_orientations.getValue();
    size_t N = omegas.size();

    // Convert stored rotation vectors → SO3 once for this cache update
    std::vector<SO3> R(N);
    for (size_t i = 0; i < N; ++i)
        R[i] = SO3::exp(SO3::TangentVector(omegas[i].x(), omegas[i].y(), omegas[i].z()));

    cached_linear_strains.assign(N, Vec3d(0.0, 0.0, 0.0));
    double max_gamma = 0.0;
    for (size_t i = 0; i < N; ++i) {
        if (i < rest_lengths.size() && rest_lengths[i] > 1e-12) {
            Vec3d dx = pos[i + 1] - pos[i];
            SO3::Vector local_dx = R[i].inverse().act(SO3::Vector(dx.x(), dx.y(), dx.z()));
            cached_linear_strains[i] = Vec3d(local_dx.x(), local_dx.y(), local_dx.z()) / rest_lengths[i] - Vec3d(1.0, 0.0, 0.0);
            max_gamma = std::max(max_gamma, cached_linear_strains[i].norm());
        }
    }

    cached_angular_strains.assign(N + 1, Vec3d(0.0, 0.0, 0.0));
    double max_omega_strain = 0.0;
    for (size_t i = 1; i < N; ++i) {
        // Angular strain (curvature + torsion) at node i, staggered convention:
        //   Omega_i = log(R_{i-1}^T * R_i) / h_dual
        // Boundary strains (i=0 and i=N) remain zero (clamped/free BCs).
        const double h_dual = (rest_lengths[i - 1] + rest_lengths[i]) * 0.5;
        if (h_dual > 1e-12) {
            SO3 relative_R = R[i - 1].inverse() * R[i];
            SO3::TangentVector omega = relative_R.log();
            cached_angular_strains[i] = Vec3d(omega.x(), omega.y(), omega.z()) / h_dual;
            max_omega_strain = std::max(max_omega_strain, cached_angular_strains[i].norm());
        }
    }

    if (printCache) {
        msg_info() << "[updateStrainsCache] max|Gamma|=" << std::scientific << std::setprecision(3)
                   << max_gamma << "  max|Omega|=" << max_omega_strain
                   << "  [both should be 0 at rest]";
        if (max_gamma > 0.5)
            msg_warning() << "[updateStrainsCache] ⚠  max|Gamma|=" << max_gamma
                          << " > 0.5 — large linear strain, possible instability";
        if (max_omega_strain > 10.0)
            msg_warning() << "[updateStrainsCache] ⚠  max|Omega|=" << max_omega_strain
                          << " rad/m — very large angular strain (beam may have flipped)";
    }
}

sofa::type::Vec3d CosseratIntrinsicState::getLinearStrain(size_t i) {
    updateStrainsCache();
    if (i < cached_linear_strains.size()) {
        return cached_linear_strains[i];
    }
    return sofa::type::Vec3d(0.0, 0.0, 0.0);
}

sofa::type::Vec3d CosseratIntrinsicState::getAngularStrain(size_t i) {
    updateStrainsCache();
    if (i < cached_angular_strains.size()) {
        return cached_angular_strains[i];
    }
    return sofa::type::Vec3d(0.0, 0.0, 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// J_r⁻¹(ω) — right Lie Jacobian inverse of SO(3)
//
// NOTE — duplication intentionnelle :
//   La même formule existe dans `sofa::component::cosserat::liegroups::SO3
//   ::computeRightJacobianInverse` (fichier src/liegroups/SO3.inl). Cette copie
//   locale est conservée pour éviter les conversions Eigen ↔ sofa::type::Mat3x3d
//   à chaque appel depuis `PainlessBeamForceField` (chemin chaud, called dans la
//   boucle de force à chaque step). Toute correction de l'une DOIT être répercutée
//   dans l'autre — les deux retournent EXACTEMENT la même fonction.
//
// Formule (Solà 2018, eq. 144) :
//   J_r⁻¹(ω) = I + ½ [ω]× + c(θ) [ω]×²
//   c(θ) = 1/θ² − (1+cos θ)/(2θ sin θ)
//   c(θ) → 1/12 quand θ → 0 (Taylor à l'ordre 3 inclus : c ≈ 1/12 + θ²/720)
// ─────────────────────────────────────────────────────────────────────────────
sofa::type::Mat3x3d CosseratIntrinsicState::getInverseLieJacobian(const sofa::type::Vec3d& omega) {
    double theta2 = omega.norm2();
    sofa::type::Mat3x3d I;
    I.identity();

    sofa::type::Mat3x3d W;
    W[0][0] = 0.0;        W[0][1] = -omega.z(); W[0][2] = omega.y();
    W[1][0] = omega.z();  W[1][1] = 0.0;        W[1][2] = -omega.x();
    W[2][0] = -omega.y(); W[2][1] = omega.x();  W[2][2] = 0.0;

    sofa::type::Mat3x3d W2 = W * W;

    if (theta2 < 1e-8) {
        // 3rd-order Taylor expansion : c(θ) ≈ 1/12 + θ²/720
        const double c = 1.0/12.0 + theta2 / 720.0;
        return I + W * 0.5 + W2 * c;
    }

    double theta = std::sqrt(theta2);
    double term2 = (1.0 / theta2) - (1.0 + std::cos(theta)) / (2.0 * theta * std::sin(theta));

    return I + W * 0.5 + W2 * term2;
}

}  // namespace sofa::component::cosserat::engine
