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

#include "Cosserat/types.h"
#include <Cosserat/forcefield/standard/HookeSeratPCSForceField.h>
#include <Cosserat/forcefield/base/HookeSeratBaseForceField.h>
#include <sofa/core/behavior/ForceField.inl>

using sofa::core::VecCoordId;
using sofa::core::behavior::MechanicalState;
using sofa::core::objectmodel::BaseContext;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;

#include <iostream>
using std::cout;
using std::endl;

#include <algorithm>
#include <ctime>

namespace sofa::component::forcefield {

	using sofa::core::behavior::BaseMechanicalState;
	using sofa::core::behavior::MultiMatrixAccessor;
	using sofa::helper::WriteAccessor;

	template<typename DataTypes>
	HookeSeratPCSForceField<DataTypes>::HookeSeratPCSForceField() :
		Inherit1(), d_variantSections(initData(&d_variantSections, false, "variantSections",
											   "In case we have variant beam sections this has to be set to true")),
		d_youngModulusList(initData(&d_youngModulusList, "youngModulusList",
									"The list of Young modulus in case we have sections with "
									"variable physical properties")),
		d_poissonRatioList(initData(&d_poissonRatioList, "poissonRatioList",
									"The list of poisson's ratio in case we have sections with "
									"variable physical properties")),
		d_useInertiaParams(initData(&d_useInertiaParams, false, "useInertiaParams",
									"If the inertia parameters are given by the user, there is no longer "
									"any need to use @d_youngModulus and @d_poissonRatio.")),
		d_GI(initData(&d_GI, "GI", "The inertia parameter, GI")),
		d_GA(initData(&d_GA, "GA", "The inertia parameter, GA")),
		d_EA(initData(&d_EA, "EA", "The inertia parameter, EA")),
		d_EI(initData(&d_EI, "EI", "The inertia parameter, EI")) {
		compute_df = true;
	}

	template<typename DataTypes>
	HookeSeratPCSForceField<DataTypes>::~HookeSeratPCSForceField() = default;

	template<typename DataTypes>
	void HookeSeratPCSForceField<DataTypes>::init() {
		Inherit1::init();
	}

	/*Cross-Section Properties Initialization: The reinit function begins by
	   recalculating the properties related to the cross-section of the beams. It
	   calculates the area moment of inertia (Iy and Iz), the polar moment of
	   inertia (J), and the cross-sectional area (A). These calculations depend on
	   the chosen cross-section shape, either circular or rectangular. T he formulas
	   used for these calculations are based on standard equations for these
	   properties.*/
	template<typename DataTypes>
	void HookeSeratPCSForceField<DataTypes>::reinit() {
		Inherit1::reinit();
		// if we are dealing with different physical properties : YM and PR
		if (!d_variantSections.getValue()) {
			if (!d_useInertiaParams.getValue()) {
				Real E = this->d_youngModulus.getValue();
				Real G = E / (2.0 * (1.0 + this->d_poissonRatio.getValue()));
				// Inertial matrix
				m_K_section(0, 0) = G * this->J;
				m_K_section(1, 1) = E * this->Iy;
				m_K_section(2, 2) = E * this->Iz;
			} else {
				msg_info("HookeSeratPCSForceField") << "Pre-calculated inertia parameters are used for the computation "
													  "of the stiffness matrix.";
				m_K_section(0, 0) = d_GI.getValue();
				m_K_section(1, 1) = d_EI.getValue();
				m_K_section(2, 2) = d_EI.getValue();
			}

		} else {
			/*If the d_variantSections flag is set to true, it implies that
			   multi-section beams are used for the simulation. In this case, the code
			   calculates and initializes a list of stiffness matrices (m_K_sectionList)
			   for each section. The properties of each section, such as Young's modulus
			   and Poisson's ratio, are specified in the d_youngModulusList and
			   d_poissonRatioList data.*/
			msg_info("HookeSeratPCSForceField") << "Multi section beam are used for the simulation!";
			m_K_sectionList.clear();

			const auto szL = this->d_length.getValue().size();
			if ((szL != d_poissonRatioList.getValue().size()) || (szL != d_youngModulusList.getValue().size())) {
				msg_error("HookeSeratPCSForceField") << "Please the size of the data length, youngModulusList and "
													   "poissonRatioList should be the same !";
				return;
			}

			/*Stiffness Matrix Initialization: Next, the code initializes the stiffness
			   matrix m_K_section based on the properties of the cross-section and the
			   material's Young's modulus (E) and Poisson's ratio. The stiffness matrix
			   is essential for computing forces and simulating beam behavior.*/
			for (size_t k = 0; k < szL; k++) {
				Matrix3 _m_K_section;
				Real E = d_youngModulusList.getValue()[k];
				Real G = E / (2.0 * (1.0 + d_poissonRatioList.getValue()[k]));

				_m_K_section(0, 0) = G * this->J;
				_m_K_section(1, 1) = E * this->Iy;
				_m_K_section(2, 2) = E * this->Iz;
				m_K_sectionList.push_back(_m_K_section);
			}
			msg_info("HookeSeratPCSForceField") << "If you plan to use a multi section beam with (different "
												  "mechanical properties) and pre-calculated inertia parameters "
												  "(GI, GA, etc.), this is not yet supported.";
		}
	}


	template<typename DataTypes>
	void HookeSeratPCSForceField<DataTypes>::addForce(const MechanicalParams *mparams, DataVecDeriv &d_f,
													 const DataVecCoord &d_x, const DataVecDeriv &d_v) {
		SOFA_UNUSED(d_v);
		SOFA_UNUSED(mparams);

		if (!this->getMState()) {
			msg_info("HookeSeratPCSForceField") << "No Mechanical State found, no force will be computed...";
			compute_df = false;
			return;
		}
		VecDeriv &f = *d_f.beginEdit();
		const VecCoord &x = d_x.getValue();
		// get the rest position (for non straight shape)
		const VecCoord &x0 = this->mstate->read(sofa::core::vec_id::read_access::restPosition)->getValue();

		f.resize(x.size());
		unsigned int sz = this->d_length.getValue().size();
		if (x.size() != sz) {
			msg_warning("HookeSeratPCSForceField")
					<< " length : " << sz << "should have the same size as x... " << x.size();
			compute_df = false;
			return;
		}

		if (!d_variantSections.getValue())
			{
			// @todo: use multithread
				for (unsigned int i = 0; i < x.size(); i++) {
					// Using the correct matrix type for the datatype
					// For Vec3Types, m_K_section should be Mat33
					Vector3 current_strain = Vector3::Map(x[i].data());
					Vector3 rest_strain = Vector3::Map(x0[i].data());
					Vector3 strain = current_strain - rest_strain;
					Vector3 _f = (m_K_section * strain) * this->d_length.getValue()[i];

					for (unsigned int j = 0; j < 3; j++)
						f[i][j] -= _f[j];
				}
	}
		else {
			// @todo: use multithread
			Vector3 current_strain, rest_strain, strain, _f;

			for (unsigned int i = 0; i < x.size(); i++) {
				current_strain = Vector3::Map(x[i].data());
				rest_strain = Vector3::Map(x0[i].data());
				strain = current_strain - rest_strain;
				_f = (m_K_sectionList[i] * strain) * this->d_length.getValue()[i];
				for (int j = 0; j < 3; j++)
					f[i][j] -= _f[j];
			}

		}

		d_f.endEdit();
	}

	template<typename DataTypes>
	void HookeSeratPCSForceField<DataTypes>::addDForce(const MechanicalParams *mparams, DataVecDeriv &d_df,
													  const DataVecDeriv &d_dx) {
		if (!compute_df)
			return;

		WriteAccessor<DataVecDeriv> df = d_df;
		ReadAccessor<DataVecDeriv> dx = d_dx;
		Real kFactor = (Real) mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
		Vector3 d_strain, _df;
		df.resize(dx.size());
		if (!d_variantSections.getValue()) {

			for (unsigned int i = 0; i < dx.size(); i++) {
				d_strain = Vector3::Map(dx[i].data());
				_df = (m_K_section * d_strain) * this->d_length.getValue()[i];
				for (unsigned int j = 0; j < 3; j++)
					df[i][j] -= _df[j];
			}
		}

		else
			for (unsigned int i = 0; i < dx.size(); i++) {
				d_strain = Vector3::Map(dx[i].data());
				_df = (m_K_sectionList[i] * d_strain) * this->d_length.getValue()[i];
				for (unsigned int j = 0; j < 3; j++)
					df[i][j] -= _df[j];
			}

	}

	template<typename DataTypes>
	double HookeSeratPCSForceField<DataTypes>::getPotentialEnergy(const MechanicalParams *mparams,
																 const DataVecCoord &d_x) const {
		SOFA_UNUSED(mparams);
		if (!this->getMState())
			return 0.0;

		const VecCoord &x = d_x.getValue();
		const VecCoord &x0 = this->mstate->read(sofa::core::vec_id::read_access::restPosition)->getValue();

		double energy = 0.0;

		if (!d_variantSections.getValue()) {
			for (unsigned int i = 0; i < x.size(); i++) {
				Vector3 strain = Vector3::Map(x[i].data()) - Vector3::Map(x0[i].data());
				// Calcul correct : 0.5 * strain^T * K * strain
				// Utilisation du produit scalaire si disponible
				energy += 0.5 * strain.dot(m_K_section * strain) * this->d_length.getValue()[i];
			}
		} else {
			for (unsigned int i = 0; i < x.size(); i++) {
				const auto &K = m_K_sectionList[i];
				const Vector3 strain = Vector3::Map(x[i].data()) - Vector3::Map(x0[i].data());

				// Utilisation du produit scalaire si disponible
				energy += 0.5 * strain.dot(K * strain) * this->d_length.getValue()[i];
			}
		}

		return energy;
	}


	template<typename DataTypes>
	void HookeSeratPCSForceField<DataTypes>::addKToMatrix(const MechanicalParams *mparams,
														 const MultiMatrixAccessor *matrix) {
		MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
		BaseMatrix *mat = mref.matrix;
		unsigned int offset = mref.offset;
		Real kFact = (Real) mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

		const VecCoord &pos = this->mstate->read(core::vec_id::read_access::position)->getValue();
		for (unsigned int n = 0; n < pos.size(); n++) {
			if (!d_variantSections.getValue())
				for (unsigned int i = 0; i < 3; i++)
					for (unsigned int j = 0; j < 3; j++)
						mat->add(offset + i + 3 * n, offset + j + 3 * n,
								 -kFact * m_K_section(i, j) * this->d_length.getValue()[n]);
			else
				for (unsigned int i = 0; i < 3; i++)
					for (unsigned int j = 0; j < 3; j++)
						mat->add(offset + i + 3 * n, offset + j + 3 * n,
								 -kFact * m_K_sectionList[n](i, j) * this->d_length.getValue()[n]);
		}
	}

	template<typename DataTypes>
	typename HookeSeratPCSForceField<DataTypes>::Real HookeSeratPCSForceField<DataTypes>::getRadius() {
		return this->d_radius.getValue();
	}

} // namespace sofa::component::forcefield