// This file defines the SE23 (extended Special Euclidean group in 3D) class,
// representing rigid body transformations with linear velocity in 3D space.

/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture * (c) 2006
 *INRIA, USTL, UJF, CNRS, MGH                     *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/\>. *
 ******************************************************************************/

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE22_H
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_H
#pragma once

#include "LieGroup.h" // Then the base class interface
#include "LieGroup.inl" // Then the base class interface
#include "SE3.h" // Then the base class interface
// #include <eigen3/Eigen/Geometry.h>

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Implementation of SE_2(3), the extended Special Euclidean group in 3D
	 *
	 * This class implements the group of rigid body transformations with linear
	 * velocity in 3D space. Elements of SE_2(3) are represented as a combination
	 * of:
	 * - An SE(3) transformation (rotation and position)
	 * - A 3D linear velocity vector
	 *
	 * The Lie algebra se_2(3) consists of vectors in ℝ⁹, where:
	 * - First three components represent linear velocity
	 * - Middle three components represent angular velocity
	 * - Last three components represent linear acceleration
	 *
	 * @tparam Scalar_ The scalar type (must be a floating-point type)
	 */
	template<typename Scalar_, int Dim_ = 9>
	class SE23 : public LieGroup<SE23<Scalar_, Dim_>, Scalar_, 9, 9, 6>
	//,public LieGroupOperations<SE23<Scalar_>>
	{
	public:
		using Base = LieGroup<SE23<Scalar_, Dim_>, Scalar_, 9, 9, 6>;
		using Scalar = typename Base::Scalar;
		using Vector = typename Base::Vector;
		using Matrix = typename Base::Matrix;
		using TangentVector = typename Base::TangentVector;
		using AdjointMatrix = typename Base::AdjointMatrix;

		using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
		using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
		using Matrix4 = Eigen::Matrix<Scalar, 4, 4>;
		static constexpr int Dim = Base::Dim;

		/**
		 * @brief Default constructor creates identity element.
		 * Initializes pose to identity and velocity to zero vector.
		 */
		SE23() : m_pose(), m_velocity(Vector3::Zero()) {}

		/**
		 * @brief Construct from pose and velocity.
		 * @param pose The SE3 pose component.
		 * @param velocity The 3D linear velocity vector.
		 */
		SE23(const SE3<Scalar> &pose, const Vector3 &velocity) : m_pose(pose), m_velocity(velocity) {}

		/**
		 * @brief Construct from rotation, position, and velocity.
		 * @param rotation The SO3 rotation component.
		 * @param position The 3D position vector.
		 * @param velocity The 3D linear velocity vector.
		 */
		SE23(const SO3<Scalar> &rotation, const Vector3 &position, const Vector3 &velocity) :
			m_pose(rotation, position), m_velocity(velocity) {}

		/**
		 * @brief Access the pose component (const version).
		 * @return A const reference to the SE3 pose component.
		 */
		const SE3<Scalar> &pose() const { return m_pose; }
		/**
		 * @brief Access the pose component (mutable version).
		 * @return A mutable reference to the SE3 pose component.
		 */
		SE3<Scalar> &pose() { return m_pose; }

		/**
		 * @brief Access the velocity component (const version).
		 * @return A const reference to the 3D linear velocity vector.
		 */
		const Vector3 &velocity() const { return m_velocity; }
		/**
		 * @brief Access the velocity component (mutable version).
		 * @return A mutable reference to the 3D linear velocity vector.
		 */
		Vector3 &velocity() { return m_velocity; }

		/**
		 * @brief Get the extended homogeneous transformation matrix.
		 * This matrix represents the SE23 element in a higher-dimensional space.
		 * @return The 5x5 extended homogeneous transformation matrix.
		 */
		Eigen::Matrix<Scalar, 5, 5> extendedMatrix() const {
			Eigen::Matrix<Scalar, 5, 5> T = Eigen::Matrix<Scalar, 5, 5>::Identity();
			T.template block<4, 4>(0, 0) = m_pose.matrix();
			T.template block<3, 1>(0, 4) = m_velocity;
			return T;
		}

		// Required CRTP methods:
		/**
		 * @brief Computes the identity element of the SE23 group.
		 * @return The identity SE23 element.
		 */
		static SE23<Scalar> computeIdentity() noexcept { return SE23(); }
		/**
		 * @brief Computes the inverse of the current SE23 element.
		 * @return The inverse SE23 element.
		 */
		SE23<Scalar> computeInverse() const {
			SE3<Scalar> inv_pose = m_pose.computeInverse();
			return SE23(inv_pose, -(inv_pose.rotation().computeAction(m_velocity)));
		}

		/**
		 * @brief Composes two SE23 elements (group multiplication).
		 * @param other The other SE23 element to compose with.
		 * @return The composed SE23 element.
		 */
		SE23<Scalar> computeCompose(const SE23 &other) const {
			SE3<Scalar> composed_pose = m_pose.compose(other.m_pose);
			Vector3 composed_velocity = m_pose.rotation().computeAction(other.m_velocity) + m_velocity;
			return SE23(composed_pose, composed_velocity);
		}

		/**
		 * @brief Computes the exponential map from the Lie algebra se_2(3) to the
		 * SE23 group.
		 * @param algebra_element The element from the Lie algebra (a 9D vector
		 * representing linear velocity, angular velocity, and linear acceleration).
		 * @return The corresponding SE23 element.
		 */
		static SE23<Scalar> computeExp(const TangentVector &algebra_element) {
			const Vector3 v = algebra_element.template segment<3>(0); // Linear velocity
			const Vector3 w = algebra_element.template segment<3>(3); // Angular velocity
			const Vector3 a = algebra_element.template segment<3>(6); // Linear acceleration

			// exp([[w^, v, a], [0, 0, 0], [0, 0, 0]]) = [[R, V v, V a], [0, 1, 0], [0, 0, 1]],
			// V(w) = I + A [w]x + B [w]x^2 being the SO(3) left Jacobian (same V for p and u). The
			// former code scaled the whole velocity by 1/theta (J a / theta with J = I + ...), so
			// the velocity column was wrong for every non-zero rotation.
			typename SE3<Scalar>::TangentVector se3_element;
			se3_element << w, v; // SE3 is angular first, linear second
			const SE3<Scalar> pose = SE3<Scalar>::computeExp(se3_element);
			Scalar A, B, C;
			detail::so3JacobianCoefficients(w.norm(), A, B, C);
			const Matrix3 W = SO3<Scalar>::computeHat(w);
			const Matrix3 V = Matrix3::Identity() + A * W + B * (W * W);
			return SE23(pose, V * a);
		}
		/**
		 * @brief Computes the logarithmic map from the SE23 group to its Lie algebra
		 * se_2(3).
		 * @return The corresponding element in the Lie algebra (a 9D vector).
		 */
		TangentVector computeLog() const {
			const typename SE3<Scalar>::TangentVector se3_part = m_pose.computeLog();
			const Vector3 w = se3_part.template head<3>(); // phi: angular velocity
			const Vector3 v = se3_part.template tail<3>(); // V^-1 p: linear velocity

			// a = V^-1 u with V^-1 = I - [w]x/2 + c [w]x^2 (half-angle coefficient, accurate up to pi)
			Scalar A, B, C;
			detail::so3JacobianCoefficients(w.norm(), A, B, C);
			const Matrix3 W = SO3<Scalar>::computeHat(w);
			const Matrix3 V_inv = Matrix3::Identity() - Scalar(0.5) * W + C * (W * W);
			TangentVector result;
			result << v, w, V_inv * m_velocity;
			return result;
		}
		/**
		 * @brief Computes the adjoint representation of the current SE23 element.
		 * @return The adjoint matrix.
		 */
		AdjointMatrix computeAdjoint() const {
			// Ad_g xi = vee(g hat(xi) g^-1), g = [[R, p, u], [0, 1, 0], [0, 0, 1]], xi = [v ; w ; a]:
			//   v' = R v + [p]x R w,  w' = R w,  a' = R a + [u]x R w
			// i.e. Ad = [[R, [p]x R, 0], [0, R, 0], [0, [u]x R, R]]. The velocity block [u]x R belongs
			// at (a, w) = (6, 3); the previous code put it at (0, 6), which maps a into v.
			AdjointMatrix Ad = AdjointMatrix::Zero();
			const Matrix3 R = m_pose.rotation().matrix();
			const Matrix3 p_hat = SO3<Scalar>::computeHat(m_pose.translation());
			const Matrix3 u_hat = SO3<Scalar>::computeHat(m_velocity);
			Ad.template block<3, 3>(0, 0) = R;
			Ad.template block<3, 3>(0, 3) = p_hat * R;
			Ad.template block<3, 3>(3, 3) = R;
			Ad.template block<3, 3>(6, 3) = u_hat * R;
			Ad.template block<3, 3>(6, 6) = R;
			return Ad;
		}
		/**
		 * @brief Checks if the current SE23 element is approximately equal to
		 * another.
		 * @param other The other SE23 element to compare with.
		 * @param eps The tolerance for approximation.
		 * @return True if the elements are approximately equal, false otherwise.
		 */
		bool computeIsApprox(const SE23 &other, const Scalar &eps = Types<Scalar>::tolerance()) const {
			// absolute tolerance scaled by max(1, |u|): Eigen's relative isApprox is false against a velocity that is 0
			const Scalar scale = std::max({Scalar(1), m_velocity.norm(), other.m_velocity.norm()});
			return m_pose.computeIsApprox(other.m_pose, eps) && (m_velocity - other.m_velocity).norm() <= eps * scale;
		}
		/**
		 * @brief Applies the group action of the current SE23 element on a
		 * point-velocity pair.
		 * @param point_vel The point-velocity pair (6D vector: 3D point, 3D
		 * velocity).
		 * @return The transformed point-velocity pair.
		 */
		typename Base::ActionVector computeAction(const typename Base::ActionVector &point_vel) const {
			Vector3 point = point_vel.template head<3>();
			Vector3 vel = point_vel.template segment<3>(3);

			// Transform position and combine velocities
			Vector3 transformed_point = m_pose.computeAction(point);
			Vector3 transformed_vel = m_pose.rotation().computeAction(vel) + m_velocity;

			typename Base::ActionVector result;
			result.resize(6);
			result << transformed_point, transformed_vel;
			return result;
		}

		/**
		 * @brief Hat operator - maps a 9D Lie algebra vector to a 5x5 matrix
		 * representation. This is a placeholder, actual implementation depends on the
		 * specific representation.
		 * @param v The 9D Lie algebra vector.
		 * @return The 5x5 matrix representation.
		 */
		static Matrix computeHat(const TangentVector &v) {
			// hat([v ; w ; a]) = [[w^, v, a], [0, 0, 0], [0, 0, 0]] (5x5), embedded in the top-left
			// corner of the 9x9 AlgebraMatrix so that matrix commutators equal the Lie bracket.
			Matrix result = Matrix::Zero();
			result.template block<3, 3>(0, 0) = SO3<Scalar>::computeHat(Vector3(v.template segment<3>(3)));
			result.template block<3, 1>(0, 3) = v.template segment<3>(0);
			result.template block<3, 1>(0, 4) = v.template segment<3>(6);
			return result;
		}

		/**
		 * @brief Vee operator - inverse of hat, maps a 5x5 matrix representation to a
		 * 9D Lie algebra vector. This is a placeholder, actual implementation depends
		 * on the specific representation.
		 * @param X The 5x5 matrix representation.
		 * @return The 9D Lie algebra vector.
		 */
		static TangentVector computeVee(const Matrix &X) {
			TangentVector result;
			result.template segment<3>(0) = X.template block<3, 1>(0, 3);
			result.template segment<3>(3) = SO3<Scalar>::computeVee(Matrix3(X.template block<3, 3>(0, 0)));
			result.template segment<3>(6) = X.template block<3, 1>(0, 4);
			return result;
		}

		/**
		 * @brief Computes the adjoint representation of a Lie algebra element for
		 * SE23. This is a placeholder, actual implementation depends on the specific
		 * representation.
		 * @param v The element of the Lie algebra in vector form.
		 * @return The adjoint matrix.
		 */
		static AdjointMatrix computeAd(const TangentVector &v) {
			// ad_xi y = [xi, y]: v' = w_x v_y - w_y v_x ..., i.e. in the [v ; w ; a] layout
			//   ad = [[w^, v^, 0], [0, w^, 0], [0, a^, w^]]
			AdjointMatrix result = AdjointMatrix::Zero();
			const Matrix3 W = SO3<Scalar>::computeHat(Vector3(v.template segment<3>(3)));
			result.template block<3, 3>(0, 0) = W;
			result.template block<3, 3>(0, 3) = SO3<Scalar>::computeHat(Vector3(v.template segment<3>(0)));
			result.template block<3, 3>(3, 3) = W;
			result.template block<3, 3>(6, 3) = SO3<Scalar>::computeHat(Vector3(v.template segment<3>(6)));
			result.template block<3, 3>(6, 6) = W;
			return result;
		}

		/**
		 * @brief Generates a random SE23 element.
		 * @tparam Generator The type of the random number generator.
		 * @param gen The random number generator.
		 * @return A random SE23 element.
		 */
		template<typename Generator>
		static SE23<Scalar> computeRandom(Generator &gen) {
			return SE23(SE3<Scalar>::computeRandom(gen), Types<Scalar>::template randomVector<3>(gen));
		}

		/**
		 * @brief Prints the SE23 element to an output stream.
		 * @param os The output stream.
		 * @return The output stream.
		 */
		std::ostream &print(std::ostream &os) const {
			os << "SE23(pose=" << m_pose << ", velocity=" << m_velocity.transpose() << ")";
			return os;
		}

		/**
		 * @brief Gets the type name of the SE23 class.
		 * @return A string view of the type name.
		 */
		static constexpr std::string_view getTypeName() { return "SE23"; }

		/**
		 * @brief Checks if the current SE23 element is valid.
		 * @return True if both pose and velocity components are valid, false
		 * otherwise.
		 */
		bool computeIsValid() const { return m_pose.computeIsValid() && m_velocity.allFinite(); }

		/**
		 * @brief Normalizes the SE23 element.
		 * Normalizes the pose component.
		 */
		void computeNormalize() { m_pose.computeNormalize(); }

	private:
		SE3<Scalar> m_pose; ///< Rigid body transformation
		Vector3 m_velocity; ///< Linear velocity
	};

} // namespace sofa::component::cosserat::liegroups

// Out-of-line free functions (operator/, Euler conversions, interpolate).
#include "SE23.inl"


// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_H
