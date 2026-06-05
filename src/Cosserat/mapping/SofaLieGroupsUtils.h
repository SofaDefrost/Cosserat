/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
 *                                                                             *
 * This file is part of the Cosserat plugin.                                   *
 ******************************************************************************/
#pragma once

#include <Eigen/Geometry>
#include <liegroups/SE3.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Quat.h>
#include <sofa/type/Vec.h>

namespace Cosserat::mapping {

	/**
	 * @brief Helpers de conversion SOFA <-> Eigen <-> SE3 (liegroups).
	 *
	 * Centralise les conventions de quaternion (SOFA: [x, y, z, w], Eigen: [w, x, y, z])
	 * et la construction d'une `liegroups::SE3<double>` à partir d'un `sofa::Rigid3d`,
	 * afin d'éliminer le copier-coller dans les *Mapping*::apply / applyJ / applyJT.
	 *
	 * Toutes ces fonctions sont `inline` et `noexcept` — aucun coût à l'appel.
	 */

	using SE3Lie    = sofa::component::cosserat::liegroups::SE3<double>;
	using SO3Lie    = SE3Lie::SO3Type;
	using Vector3d  = SE3Lie::Vector3;

	/**
	 * @brief Convertit un quaternion SOFA `[x, y, z, w]` en `Eigen::Quaterniond` `[w, x, y, z]`.
	 *
	 * @note `Eigen::Quaternion` est construit dans l'ordre `(w, x, y, z)`.
	 *       SOFA `Quat<T>::operator[]` renvoie `[x, y, z, w]`.
	 */
	template<class Scalar>
	[[nodiscard]] inline Eigen::Quaternion<double>
	sofaQuatToEigen(const sofa::type::Quat<Scalar>& q) noexcept {
		return Eigen::Quaternion<double>(
			static_cast<double>(q[3]),   // w
			static_cast<double>(q[0]),   // x
			static_cast<double>(q[1]),   // y
			static_cast<double>(q[2]));  // z
	}

	/**
	 * @brief Convertit un `Eigen::Quaterniond` en quaternion SOFA `[x, y, z, w]`.
	 */
	[[nodiscard]] inline sofa::type::Quat<SReal>
	eigenToSofaQuat(const Eigen::Quaternion<double>& q) noexcept {
		return sofa::type::Quat<SReal>(
			static_cast<SReal>(q.x()),
			static_cast<SReal>(q.y()),
			static_cast<SReal>(q.z()),
			static_cast<SReal>(q.w()));
	}

	/**
	 * @brief Convertit un `sofa::type::Vec3<Scalar>` (centre d'un Rigid3) en `Vector3d` Eigen.
	 */
	template<class Scalar>
	[[nodiscard]] inline Vector3d
	sofaVec3ToEigen(const sofa::type::Vec<3, Scalar>& v) noexcept {
		return Vector3d(static_cast<double>(v[0]),
						static_cast<double>(v[1]),
						static_cast<double>(v[2]));
	}

	/**
	 * @brief Construit un `liegroups::SE3<double>` depuis un `sofa::Rigid3<Scalar>::Coord` (centre + quat).
	 *
	 * Centralise la conversion répétée :
	 *   - centre SOFA -> `Vector3d`
	 *   - quaternion SOFA `[x,y,z,w]` -> `Eigen::Quaterniond` -> `SO3`
	 */
	template<class RigidCoord>
	[[nodiscard]] inline SE3Lie
	rigidCoordToSE3(const RigidCoord& coord) noexcept {
		return SE3Lie(SO3Lie(sofaQuatToEigen(coord.getOrientation()).toRotationMatrix()),
					  sofaVec3ToEigen(coord.getCenter()));
	}

	/**
	 * @brief Inverse de `rigidCoordToSE3` : remplit un Rigid3 SOFA depuis un SE3 liegroups.
	 */
	template<class RigidCoord>
	inline void se3ToRigidCoord(const SE3Lie& g, RigidCoord& out) noexcept {
		const auto& t = g.translation();
		out.getCenter() = sofa::type::Vec<3, SReal>(
			static_cast<SReal>(t[0]),
			static_cast<SReal>(t[1]),
			static_cast<SReal>(t[2]));
		out.getOrientation() = eigenToSofaQuat(Eigen::Quaternion<double>(g.rotation().matrix()));
	}

} // namespace Cosserat::mapping
