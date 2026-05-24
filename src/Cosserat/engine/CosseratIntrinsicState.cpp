#include <Cosserat/engine/CosseratIntrinsicState.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/logging/Messaging.h>

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

void CosseratIntrinsicState::init() {
    BaseMechanicalState::init();

    const auto& pos = d_positions.getValue();
    const auto& R = d_orientations.getValue();
    size_t nbSections = R.size();

    rest_lengths.clear();
    rest_lengths.reserve(nbSections);
    
    for (size_t i = 0; i < nbSections; ++i) {
        if (i + 1 < pos.size()) {
            double L = (pos[i+1] - pos[i]).norm();
            rest_lengths.push_back(L);
        } else {
            msg_warning() << "Node indices out of bounds while computing rest lengths!";
            rest_lengths.push_back(1.0);
        }
    }
}

void CosseratIntrinsicState::updateStrainsCache() {
    bool dirty = false;
    if (d_positions.getCounter() != last_positions_counter) {
        dirty = true;
        last_positions_counter = d_positions.getCounter();
    }
    if (d_orientations.getCounter() != last_orientations_counter) {
        dirty = true;
        last_orientations_counter = d_orientations.getCounter();
    }

    if (!dirty) return;

    const auto& pos = d_positions.getValue();
    const auto& R = d_orientations.getValue();
    size_t N = R.size();

    cached_linear_strains.assign(N, Vec3d(0.0, 0.0, 0.0));
    for (size_t i = 0; i < N; ++i) {
        if (i < rest_lengths.size() && rest_lengths[i] > 1e-12) {
            Vec3d dx = pos[i + 1] - pos[i];
            SO3::Vector local_dx = R[i].inverse().act(SO3::Vector(dx.x(), dx.y(), dx.z()));
            cached_linear_strains[i] = Vec3d(local_dx.x(), local_dx.y(), local_dx.z()) / rest_lengths[i] - Vec3d(1.0, 0.0, 0.0);
        }
    }

    cached_angular_strains.assign(N + 1, Vec3d(0.0, 0.0, 0.0));
    for (size_t i = 1; i < N; ++i) {
        // Angular strain (curvature + torsion) at node i, staggered convention:
        //   Omega_i = log(R_{i-1}^T * R_i) / h_dual
        // where h_dual is the length of the dual edge connecting the midpoints of
        // segments i-1 and i:
        //   h_dual = (h_{i-1} + h_i) / 2
        // For a uniform beam all h_i are equal so h_dual = h; the general case
        // requires averaging the two adjacent segment lengths.
        // Boundary strains (i=0 and i=N) remain zero (clamped/free BCs applied
        // at the ForceField level, not here).
        const double h_dual = (rest_lengths[i - 1] + rest_lengths[i]) * 0.5;
        if (h_dual > 1e-12) {
            SO3 relative_R = R[i - 1].inverse() * R[i];
            SO3::TangentVector omega = relative_R.log();
            cached_angular_strains[i] = Vec3d(omega.x(), omega.y(), omega.z()) / h_dual;
        }
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

sofa::type::Mat3x3d CosseratIntrinsicState::getInverseLieJacobian(const sofa::type::Vec3d& omega) {
    double theta2 = omega.norm2();
    sofa::type::Mat3x3d I = sofa::type::Mat3x3d::identity();
    
    sofa::type::Mat3x3d W;
    W[0][0] = 0.0;        W[0][1] = -omega.z(); W[0][2] = omega.y();
    W[1][0] = omega.z();  W[1][1] = 0.0;        W[1][2] = -omega.x();
    W[2][0] = -omega.y(); W[2][1] = omega.x();  W[2][2] = 0.0;
    
    sofa::type::Mat3x3d W2 = W * W;

    if (theta2 < 1e-12) {
        // Taylor expansion for small theta
        return I + W * 0.5 + W2 * (1.0 / 12.0);
    }

    double theta = std::sqrt(theta2);
    double term2 = (1.0 / theta2) - (1.0 + std::cos(theta)) / (2.0 * theta * std::sin(theta));

    return I + W * 0.5 + W2 * term2;
}

}  // namespace sofa::component::cosserat::engine
