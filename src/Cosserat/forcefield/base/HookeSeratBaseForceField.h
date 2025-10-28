#pragma once

#include <Cosserat/config.h>

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/linearalgebra/CompressedRowSparseMatrix.h>

#include <liegroups/Types.h>

namespace sofa::component::forcefield {

	using sofa::core::MechanicalParams;
	using sofa::core::behavior::ForceField;
	using sofa::core::behavior::MultiMatrixAccessor;
	using sofa::helper::OptionsGroup;
	using sofa::linearalgebra::BaseMatrix;
	using sofa::linearalgebra::CompressedRowSparseMatrix;
	using sofa::type::Mat;
	using sofa::type::Vec;
	using sofa::type::vector;

	/**
	 * @brief Base class for beam force field implementations
	 *
	 * This class provides the common functionality for beam-based force fields,
	 * including cross-section property calculations and material parameters.
	 *
	 * @tparam DataTypes The data types used for coordinates and derivatives
	 */
	template<typename DataTypes>
	class HookeSeratBaseForceField : public ForceField<DataTypes> {
	public:
		SOFA_CLASS(SOFA_TEMPLATE(HookeSeratBaseForceField, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));

		// Type definitions
		using Real = typename DataTypes::Real;
		using VecCoord = typename DataTypes::VecCoord;
		using VecDeriv = typename DataTypes::VecDeriv;
		using Coord = typename DataTypes::Coord;
		using Deriv = typename DataTypes::Deriv;

		using DataVecCoord = Data<VecCoord>;
		using DataVecDeriv = Data<VecDeriv>;

		// Encapsulate Lie group types
		struct LieTypes {
			using Types = sofa::component::cosserat::liegroups::Types<double>;
			using Vector2 = typename Types::Vector2;
			using Vector3 = typename Types::Vector3;
			using Vector4 = typename Types::Vector4;
			using Vector6 = typename Types::Vector6;
			using Matrix3 = typename Types::Matrix3;
			using Matrix6 = typename Types::Matrix6;

		};

		using Vector2 = typename LieTypes::Vector2;
		using Vector3 = typename LieTypes::Vector3;
		using Vector4 = typename LieTypes::Vector4;
		using Vector6 = typename LieTypes::Vector6;
		using Matrix3 = typename LieTypes::Matrix3;
		using Matrix6 = typename LieTypes::Matrix6;

		/**
		 * @brief Enumeration for cross-section shapes
		 */
		enum class CrossSectionShape { CIRCULAR = 0, RECTANGULAR = 1, HOLLOW_CIRCULAR = 2 };

	public:
		HookeSeratBaseForceField();
		virtual ~HookeSeratBaseForceField() = default;

		////////////////////////// Inherited from BaseObject /////////////////////////
		void init() override;
		void reinit() override;
		///////////////////////////////////////////////////////////////////////////

		////////////////////////// Getters for derived classes ////////////////////
		/**
		 * @brief Get the cross-section area
		 * @return The cross-section area
		 */
		Real getCrossSectionArea() const { return m_crossSectionArea; }

		/**
		 * @brief Get the torsional constant J
		 * @return The torsional constant
		 */
		Real getTorsionalConstant() const { return J; }

		/**
		 * @brief Get the second moment of area about y-axis
		 * @return Second moment of area Iy
		 */
		Real getSecondMomentY() const { return Iy; }

		/**
		 * @brief Get the second moment of area about z-axis
		 * @return Second moment of area Iz
		 */
		Real getSecondMomentZ() const { return Iz; }

		/**
		 * @brief Get the Young's modulus
		 * @return Young's modulus
		 */
		Real getYoungModulus() const { return d_youngModulus.getValue(); }

		/**
		 * @brief Get the Poisson's ratio
		 * @return Poisson's ratio
		 */
		Real getPoissonRatio() const { return d_poissonRatio.getValue(); }

		/**
		 * @brief Get the shear modulus (calculated from E and ν)
		 * @return Shear modulus G
		 */
		Real getShearModulus() const { return d_youngModulus.getValue() / (2.0 * (1.0 + d_poissonRatio.getValue())); }

		/**
		 * @brief Get the current cross-section shape
		 * @return Cross-section shape enumeration
		 */
		CrossSectionShape getCrossSectionShape() const {
			return static_cast<CrossSectionShape>(d_crossSectionShape.getValue().getSelectedId());
		}

		/**
		 * @brief Check if the beam configuration is valid
		 * @return True if valid, false otherwise
		 */
		bool isValidConfiguration() const;
		///////////////////////////////////////////////////////////////////////////

		// Convert Coord to Vector6 type
		Vector6 convert_to_local_type(const Coord & x) const {
			Vector6 _x;
			for (unsigned int i = 0; i < 6; i++) {
				_x(i) = x[i];
			}
			return _x;
		}

	protected:
		using ForceField<DataTypes>::mstate;
		////////////////////////// Data members //////////////////////////////////
		/// Cross-section shape selector
		Data<helper::OptionsGroup> d_crossSectionShape;

		/// Material properties
		Data<Real> d_youngModulus;
		Data<Real> d_poissonRatio;

		/// Geometric properties
		Data<type::vector<Real>> d_length;

		/// Circular Cross Section parameters
		Data<Real> d_radius;
		Data<Real> d_innerRadius; ///< For hollow circular sections

		/// Rectangular Cross Section parameters
		Data<Real> d_lengthY; ///< Width in Y direction
		Data<Real> d_lengthZ; ///< Height in Z direction
		///////////////////////////////////////////////////////////////////////////

		////////////////////////// Computed properties ///////////////////////////
		/// Cross-section area
		Real m_crossSectionArea;

		/// Torsional constant (polar moment of inertia)
		Real J;

		/// Second moments of area
		Real Iy; ///< Second moment of area about y-axis
		Real Iz; ///< Second moment of area about z-axis
		///////////////////////////////////////////////////////////////////////////

		////////////////////////// Protected methods /////////////////////////////
		/**
		 * @brief Compute cross-section properties based on shape and dimensions
		 */
		virtual void computeCrossSectionProperties();

		/**
		 * @brief Compute properties for circular cross-section
		 */
		void computeCircularProperties();

		/**
		 * @brief Compute properties for rectangular cross-section
		 */
		void computeRectangularProperties();

		/**
		 * @brief Compute properties for hollow circular cross-section
		 */
		void computeHollowCircularProperties();

		/**
		 * @brief Validate input parameters
		 * @return True if all parameters are valid
		 */
		bool validateParameters() const;

		/**
		 * @brief Print debug information about computed properties
		 */
		void printDebugInfo() const;
		///////////////////////////////////////////////////////////////////////////


	private:
		////////////////////////// Private helper methods ////////////////////////
		/**
		 * @brief Initialize cross-section shape options
		 */
		void initializeCrossSectionOptions();

		/**
		 * @brief Check if geometric parameters are consistent
		 * @return True if consistent, false otherwise
		 */
		bool areGeometricParametersConsistent() const;
		///////////////////////////////////////////////////////////////////////////

		////////////////////////// Inherited attributes (for lookup) /////////////
		/// Bring inherited attributes into current lookup context
		using ForceField<DataTypes>::getContext;
		using ForceField<DataTypes>::f_printLog;

		///////////////////////////////////////////////////////////////////////////
	};

////////////////////////// External template declarations ////////////////////
#if !defined(SOFA_COSSERAT_CPP_HookeSeratBaseForceField)
	extern template class SOFA_COSSERAT_API HookeSeratBaseForceField<defaulttype::Vec3Types>;
	extern template class SOFA_COSSERAT_API HookeSeratBaseForceField<defaulttype::Vec6Types>;
#endif

} // namespace sofa::component::forcefield
