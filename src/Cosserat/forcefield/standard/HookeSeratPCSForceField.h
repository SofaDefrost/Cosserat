/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture                          *
 *                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
 *                                                                             *
 * This library is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this library; if not, write to the Free Software Foundation,     *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
 *******************************************************************************
 *                           Plugin Cosserat    v1.0                           *
 *				                                              *
 * This plugin is also distributed under the GNU LGPL (Lesser General          *
 * Public License) license with the same conditions than SOFA.                 *
 *                                                                             *
 * Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
 *               Ecole Centrale de Lille)                                      *
 *                                                                             *
 * Contact information: https://project.inria.fr/softrobot/contact/            *
 *                                                                             *
 ******************************************************************************/
#pragma once
#include <Cosserat/forcefield/base/HookeSeratBaseForceField.h>

namespace sofa::component::forcefield {

	/**
	 * This component is used to compute the Hooke's law on a beam computed on strain / stress
	 * Only bending and torsion strain / stress are considered here for Vec3Types
	 * Full 6DOF strain/stress for Vec6Types
	 */
	template<typename DataTypes>
	class HookeSeratPCSForceField : public HookeSeratBaseForceField<DataTypes> {
	public:
		SOFA_CLASS(SOFA_TEMPLATE(HookeSeratPCSForceField, DataTypes), SOFA_TEMPLATE(HookeSeratBaseForceField, DataTypes));

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

	public:
		HookeSeratPCSForceField();
		virtual ~HookeSeratPCSForceField();

		////////////////////////// Inherited from BaseObject /////////////////////////
		void init() override;
		void reinit() override;
		///////////////////////////////////////////////////////////////////////////

		////////////////////////// Inherited from ForceField /////////////////////////
		void addForce(const MechanicalParams *mparams, DataVecDeriv &f, const DataVecCoord &x,
					  const DataVecDeriv &v) override;

		void addDForce(const MechanicalParams *mparams, DataVecDeriv &df, const DataVecDeriv &dx) override;


		void addKToMatrix(const MechanicalParams *mparams, const MultiMatrixAccessor *matrix) override;

		double getPotentialEnergy(const MechanicalParams *mparams, const DataVecCoord &x) const override;
		////////////////////////////////////////////////////////////////////////////

		Real getRadius();

	protected:
		// In case we have a beam with different properties per section
		Data<bool> d_variantSections; /// bool to identify different beams sections
		Data<type::vector<Real>> d_youngModulusList; /// youngModulus
		Data<type::vector<Real>> d_poissonRatioList; /// poissonRatio
		/// If the inertia parameters are given by the user, there is no longer any need to use YG.
		Data<bool> d_useInertiaParams;
		Data<Real> d_GI;
		Data<Real> d_GA;
		Data<Real> d_EA;
		Data<Real> d_EI;
		Data<Real> d_EIy;
		Data<Real> d_EIz;

		bool compute_df;

		// The stiffness matrix for the beam section
		Mat33 m_K_section;
		type::vector<Mat33> m_K_sectionList;

		// The stiffness matrix for the beam section in 6x6 format
		Mat66 m_K_section66;
		type::vector<Mat66> m_k_section66List;

	private:
		////////////////////////// Inherited attributes ////////////////////////////
		/// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
		/// Bring inherited attributes and function in the current lookup context.
		/// otherwise any access to the base::attribute would require
		/// the "this->" approach.
		using ForceField<DataTypes>::getContext;
		using ForceField<DataTypes>::f_printLog;
		////////////////////////////////////////////////////////////////////////////
	};

	// Explicit declaration of this spécialisation
	template<>
	void HookeSeratPCSForceField<sofa::defaulttype::Vec6Types>::reinit();
	template<>
	void HookeSeratPCSForceField<sofa::defaulttype::Vec6Types>::addForce(const MechanicalParams *mparams,
																		 DataVecDeriv &f, const DataVecCoord &x,
																		 const DataVecDeriv &v);
	template<>
	void HookeSeratPCSForceField<sofa::defaulttype::Vec6Types>::addDForce(const MechanicalParams *mparams,
																		  DataVecDeriv &df, const DataVecDeriv &dx);
	template<>
	void HookeSeratPCSForceField<sofa::defaulttype::Vec6Types>::addKToMatrix(const MechanicalParams *mparams,
																			 const MultiMatrixAccessor *matrix);
	template<>
	double HookeSeratPCSForceField<sofa::defaulttype::Vec6Types>::getPotentialEnergy(const MechanicalParams *mparams,
																					 const DataVecCoord &x) const;

#if !defined(SOFA_COSSERAT_CPP_HookeSeratPCSForceField)
	extern template class SOFA_COSSERAT_API HookeSeratPCSForceField<defaulttype::Vec3Types>;
	extern template class SOFA_COSSERAT_API HookeSeratPCSForceField<defaulttype::Vec6Types>;
#endif

} // namespace sofa::component::forcefield
