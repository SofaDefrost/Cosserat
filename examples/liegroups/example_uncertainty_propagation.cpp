// Example demonstrating uncertainty propagation on Lie groups
// This example shows how to use the GaussianOnManifold class for
// state estimation applications.

#include <Eigen/Dense>
#include <iostream>
#include <vector>

// Forward declarations for Lie group types
namespace sofa::component::cosserat::liegroups {

	template<typename Scalar>
	class SO3 {
	public:
		using TangentVector = Eigen::Matrix<Scalar, 3, 1>;
		using AdjointMatrix = Eigen::Matrix<Scalar, 3, 3>;

		SO3() = default;
		static SO3 exp(const TangentVector &omega);
		TangentVector log() const;
		AdjointMatrix adjoint() const;
		SO3 inverse() const;
		SO3 plus(const TangentVector &tau) const { return compose(exp(tau)); }
		TangentVector minus(const SO3 &other) const { return inverse().compose(other).log(); }
		SO3 compose(const SO3 &other) const;
		static AdjointMatrix rightJacobian(const TangentVector &tau);
	};

	template<typename Scalar>
	class SE3 {
	public:
		using TangentVector = Eigen::Matrix<Scalar, 6, 1>;
		using AdjointMatrix = Eigen::Matrix<Scalar, 6, 6>;

		SE3() = default;
		static SE3 exp(const TangentVector &xi);
		TangentVector log() const;
		AdjointMatrix adjoint() const;
		SE3 inverse() const;
		SE3 plus(const TangentVector &tau) const { return compose(exp(tau)); }
		TangentVector minus(const SE3 &other) const { return inverse().compose(other).log(); }
		SE3 compose(const SE3 &other) const;
		static AdjointMatrix rightJacobian(const TangentVector &tau);
	};

	// Simplified Gaussian distribution on Lie groups
	template<typename LieGroupType>
	class GaussianOnManifold {
	public:
		using Scalar = typename LieGroupType::Scalar;
		using TangentVector = typename LieGroupType::TangentVector;
		using Covariance = typename LieGroupType::AdjointMatrix;

		GaussianOnManifold(const LieGroupType &mean, const Covariance &covariance) :
			m_mean(mean), m_covariance(covariance) {}

		const LieGroupType &mean() const { return m_mean; }
		const Covariance &covariance() const { return m_covariance; }

		// Transform covariance between frames
		Covariance toGlobalFrame() const {
			auto Ad = m_mean.adjoint();
			return Ad * m_covariance * Ad.transpose();
		}

		// Propagate through composition
		GaussianOnManifold<LieGroupType> composeWith(const TangentVector &delta, const Covariance &delta_cov) const {

			auto new_mean = m_mean.plus(delta);
			auto F = LieGroupType::exp(delta).inverse().adjoint();
			auto G = LieGroupType::rightJacobian(delta);
			auto new_cov = F * m_covariance * F.transpose() + G * delta_cov * G.transpose();

			return GaussianOnManifold<LieGroupType>(new_mean, new_cov);
		}

	private:
		LieGroupType m_mean;
		Covariance m_covariance;
	};

} // namespace sofa::component::cosserat::liegroups

using namespace sofa::component::cosserat::liegroups;

int main() {
	std::cout << "Lie Groups Uncertainty Propagation Example" << std::endl;
	std::cout << "==========================================" << std::endl;

	// Example 1: SO(3) uncertainty propagation
	std::cout << "\n1. SO(3) Rotation Uncertainty Propagation" << std::endl;
	std::cout << "----------------------------------------" << std::endl;

	// Create initial rotation with uncertainty
	SO3<double> initial_rotation; // Identity rotation
	Eigen::Matrix3d initial_cov = Eigen::Matrix3d::Identity() * 0.01; // Small uncertainty
	GaussianOnManifold<SO3<double>> rotation_dist(initial_rotation, initial_cov);

	std::cout << "Initial rotation covariance trace: " << initial_cov.trace() << std::endl;

	// Apply a rotation increment with its own uncertainty
	Eigen::Vector3d rotation_increment(0.1, 0.05, 0.02); // Small rotation
	Eigen::Matrix3d increment_cov = Eigen::Matrix3d::Identity() * 0.005; // Increment uncertainty

	auto propagated_rotation = rotation_dist.composeWith(rotation_increment, increment_cov);

	std::cout << "Propagated rotation covariance trace: " << propagated_rotation.covariance().trace() << std::endl;

	// Example 2: SE(3) pose uncertainty propagation
	std::cout << "\n2. SE(3) Pose Uncertainty Propagation" << std::endl;
	std::cout << "-------------------------------------" << std::endl;

	// Create initial pose with uncertainty
	SE3<double> initial_pose; // Identity pose
	Eigen::Matrix<double, 6, 6> pose_cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.01;
	GaussianOnManifold<SE3<double>> pose_dist(initial_pose, pose_cov);

	std::cout << "Initial pose covariance trace: " << pose_cov.trace() << std::endl;

	// Apply a pose increment (translation + rotation)
	Eigen::Matrix<double, 6, 1> pose_increment;
	pose_increment << 0.1, 0.05, 0.02, 0.01, 0.005, 0.003; // [vx, vy, vz, ωx, ωy, ωz]
	Eigen::Matrix<double, 6, 6> pose_increment_cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.002;

	auto propagated_pose = pose_dist.composeWith(pose_increment, pose_increment_cov);

	std::cout << "Propagated pose covariance trace: " << propagated_pose.covariance().trace() << std::endl;

	// Example 3: Frame transformations
	std::cout << "\n3. Covariance Frame Transformations" << std::endl;
	std::cout << "-----------------------------------" << std::endl;

	// Transform covariance to global frame
	auto global_cov = rotation_dist.toGlobalFrame();
	std::cout << "Local frame covariance trace: " << rotation_dist.covariance().trace() << std::endl;
	std::cout << "Global frame covariance trace: " << global_cov.trace() << std::endl;

	std::cout << "\nExample completed successfully!" << std::endl;
	std::cout << "This demonstrates the basic concepts of uncertainty propagation on Lie groups." << std::endl;
	std::cout << "In a real implementation, you would use the actual Lie group classes with proper Jacobians."
			  << std::endl;

	return 0;
}
