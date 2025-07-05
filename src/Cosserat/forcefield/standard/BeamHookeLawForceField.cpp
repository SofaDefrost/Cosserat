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
#define SOFA_COSSERAT_CPP_BeamHookeLawForceField

#include <Cosserat/forcefield/standard/BeamHookeLawForceField.h>
#include <Cosserat/forcefield/standard/BeamHookeLawForceField.inl>

#include <sofa/core/ObjectFactory.h>

using namespace sofa::defaulttype;

namespace sofa::component::forcefield {

// // Implementation for Vec6Types
// template class BeamHookeLawForceField<defaulttype::Vec6Types>;
//
// // Define specializations before the template instantiation
// 	template <>
// void BeamHookeLawForceField<defaulttype::Vec6Types>::reinit() {
// 		Inherit1::reinit();
// 		if (d_useInertiaParams.getValue()) {
// 			msg_info() << "Pre-calculated inertia parameters are used for the computation of the stiffness matrix.";
// 			m_K_section66[0][0] = d_GI.getValue();
// 			m_K_section66[1][1] = d_EIy.getValue();
// 			m_K_section66[2][2] = d_EIz.getValue();
// 			m_K_section66[3][3] = d_EA.getValue();
// 			m_K_section66[4][4] = d_GA.getValue();
// 			m_K_section66[5][5] = d_GA.getValue();
// 		} else {
// 			// Pour du vec 6, on a  _m_K =diag([G*J E*Iy E*Iz E*A G*As G*As]); % stifness matrix
// 			Real E = this->d_youngModulus.getValue();
// 			Real G = E / (2.0 * (1.0 + this->d_poissonRatio.getValue()));
// 			// Translational Stiffness (X, Y, Z directions):
// 			//  Gs * J(i): Shearing modulus times the second moment of the area (torsional stiffness). E * Js(i):
// 			//  Young's modulus times the second moment of the area (bending stiffness).
// 			m_K_section66[0][0] = G * this->J;
// 			m_K_section66[1][1] = E * this->Iy;
// 			m_K_section66[2][2] = E * this->Iz;
// 			// Rotational Stiffness (X, Y, Z directions):
// 			// E * A: Young's modulus times the cross-sectional area (axial stiffness).
// 			// Gs * A: Shearing modulus times the cross-sectional area (shearing stiffness).
// 			m_K_section66[3][3] = E * this->m_crossSectionArea;
// 			m_K_section66[4][4] = G * this->m_crossSectionArea;
// 			m_K_section66[5][5] = G * this->m_crossSectionArea;
// 		}
// 	}
//
// 	template <>
// void BeamHookeLawForceField<defaulttype::Vec6Types>::addForce(const MechanicalParams *mparams, DataVecDeriv &d_f,
// 																  const DataVecCoord &d_x, const DataVecDeriv &d_v) {
// 		SOFA_UNUSED(d_v);
// 		SOFA_UNUSED(mparams);
//
// 		if (!this->getMState()) {
// 			msg_info() << "No Mechanical State found, no force will be computed..." << "\n";
// 			compute_df = false;
// 			return;
// 		}
// 		VecDeriv &f = *d_f.beginEdit();
// 		const VecCoord &x = d_x.getValue();
// 		// get the rest position (for non straight shape)
// 		const VecCoord &x0 = this->mstate->read(sofa::core::vec_id::read_access::restPosition)->getValue();
//
// 		f.resize(x.size());
// 		unsigned int sz = this->d_length.getValue().size();
// 		if (x.size() != sz) {
// 			msg_warning() << " length : " << sz << "should have the same size as x... " << x.size() << "\n";
// 			compute_df = false;
// 			return;
// 		}
// 		for (unsigned int i = 0; i < x.size(); i++) {
// 			f[i] -= (m_K_section66 * (x[i] - x0[i])) * this->d_length.getValue()[i];
// 		}
//
// 		d_f.endEdit();
// 	}
//
// 	template <>
// void BeamHookeLawForceField<defaulttype::Vec6Types>::addDForce(const MechanicalParams *mparams, DataVecDeriv &d_df,
// 																   const DataVecDeriv &d_dx) {
// 		if (!compute_df)
// 			return;
//
// 		WriteAccessor<DataVecDeriv> df = d_df;
// 		ReadAccessor<DataVecDeriv> dx = d_dx;
// 		const Real kFactor =
// 				static_cast<Real>(mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue()));
//
// 		df.resize(dx.size());
// 		for (unsigned int i = 0; i < dx.size(); i++) {
// 			df[i] -= (m_K_section66 * dx[i]) * kFactor * this->d_length.getValue()[i];
// 		}
// 	}

// 	template <>
// void BeamHookeLawForceField<defaulttype::Vec6Types>::addKToMatrix(const MechanicalParams *mparams,
// 																	  const MultiMatrixAccessor *matrix) {
// 		MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
// 		BaseMatrix *mat = mref.matrix;
// 		unsigned int offset = mref.offset;
// 		const Real kFact = (Real) mparams->kFactorIncludingRayleighDamping(-this->rayleighStiffness.getValue());
// 		// static_cast<Real>(mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue()));
//
// 		const VecCoord &pos = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
// 		constexpr int VEC_SIZE = 6;
// 		for (unsigned int n = 0; n < pos.size(); n++) {
// 			for (unsigned int i = 0; i < VEC_SIZE; i++) {
// 				for (unsigned int j = 0; j < VEC_SIZE; j++) {
// 					mat->add(offset + i + VEC_SIZE * n, offset + j + VEC_SIZE * n,
// 							 -kFact * m_K_section66[i][j] * this->d_length.getValue()[n]);
// 				}
// 			}
// 		}
// 	}

	// template<>
	// double BeamHookeLawForceField<defaulttype::Vec6Types>::getPotentialEnergy(const MechanicalParams *mparams,
	// 																		  const DataVecCoord &d_x) const {
	// 	SOFA_UNUSED(mparams);
	// 	if (!this->getMState())
	// 		return 0.0;
	//
	// 	const VecCoord &x = d_x.getValue();
	// 	const VecCoord &x0 = this->mstate->read(sofa::core::vec_id::read_access::restPosition)->getValue();
	//
	// 	double energy = 0.0;
	// 	for (unsigned int i = 0; i < x.size(); i++) {
	// 		const auto &K = m_K_section66;
	// 		const auto &strain = x[i] - x0[i];
	// 		energy += 0.5 * (strain * (K * strain)) * this->d_length.getValue()[i];
	// 	}
	// 	return energy;
	// }

// Explicit template Instantiation for Vec3Types
	template class SOFA_COSSERAT_API BeamHookeLawForceField<Vec3Types>;
//template class SOFA_COSSERAT_API BeamHookeLawForceField< sofa::defaulttype::Vec6Types>;

} // namespace sofa::component::forcefield

using namespace sofa::defaulttype;

namespace Cosserat {

	void registerBeamHookeLawForceField(sofa::core::ObjectFactory *factory) {
		factory->registerObjects(
				sofa::core::ObjectRegistrationData(
						"This component is used to compute internal stress for torsion (along x) and bending (y and z)")
						.add<sofa::component::forcefield::BeamHookeLawForceField<sofa::defaulttype::Vec3Types>>(true));
						//.add<sofa::component::forcefield::BeamHookeLawForceField<sofa::defaulttype::Vec6Types>>());
	}

} // namespace Cosserat
