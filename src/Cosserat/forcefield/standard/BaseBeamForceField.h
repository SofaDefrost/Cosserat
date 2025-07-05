#pragma once

#include <Cosserat/config.h>

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/linearalgebra/CompressedRowSparseMatrix.h>

// Include the liegroups Types.h for Matrix3
#include "../../../liegroups/Types.h"

namespace sofa::component::forcefield {

	using sofa::core::MechanicalParams;
	using sofa::core::behavior::ForceField;
	using sofa::core::behavior::MultiMatrixAccessor;
	using sofa::linearalgebra::BaseMatrix;
	using sofa::linearalgebra::CompressedRowSparseMatrix;
using sofa::type::Mat;
using sofa::type::Vec;
using sofa::type::vector;

using sofa::helper::OptionsGroup;

// Using types from liegroups
using sofa::component::cosserat::liegroups::Typesd;

	template<typename DataTypes>
	class BaseBeamForceField : public ForceField<DataTypes> {
	public:
		SOFA_CLASS(SOFA_TEMPLATE(BaseBeamForceField, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));

		typedef typename DataTypes::Real Real;
		typedef typename DataTypes::VecCoord VecCoord;
		typedef typename DataTypes::VecDeriv VecDeriv;
		typedef typename DataTypes::Coord Coord;
		typedef typename DataTypes::Deriv Deriv;

		typedef Data<VecCoord> DataVecCoord;
		typedef Data<VecDeriv> DataVecDeriv;

	typedef Vec<3, Real> Vec3;
	typedef Mat<3, 3, Real> Mat33;
	typedef Mat<6, 6, Real> Mat66;
	// Use Matrix3 from liegroups Types
	typedef Typesd::Matrix3 Matrix3;

	public:
		BaseBeamForceField();
		virtual ~BaseBeamForceField() = default;

		void init() override;
		void reinit() override;

	protected:
		Data<helper::OptionsGroup> d_crossSectionShape;
		Data<Real> d_youngModulus;
		Data<Real> d_poissonRatio;
		Data<type::vector<Real>> d_length;

		/// Circular Cross Section
		Data<Real> d_radius;
		Data<Real> d_innerRadius;
		/// Rectangular Cross Section
		Data<Real> d_lengthY;
		Data<Real> d_lengthZ;

	protected:
		/// Cross-section area
		Real m_crossSectionArea;
		Real J, Iy, Iz;
	};

} // namespace sofa::component::forcefield
