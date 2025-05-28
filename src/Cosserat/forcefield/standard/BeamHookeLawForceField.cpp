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
#include <Cosserat/forcefield/BeamHookeLawForceField.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa::component::forcefield
{

template<>
void BeamHookeLawForceField<defaulttype::Vec6Types>::addForce(
    const MechanicalParams* mparams,
    DataVecDeriv& d_f,
    const DataVecCoord& d_x,
    const DataVecDeriv& d_v)
{
    SOFA_UNUSED(d_v);
    SOFA_UNUSED(mparams);

    if (!this->m_initialized) {
        msg_error() << "Force field not properly initialized";
        return;
    }

    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    const VecCoord& x0 = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
    const vector<Real>& lengths = this->d_length.getValue();

    f.resize(x.size());

    // Compute forces using LieGroups
    for (size_t i = 0; i < x.size(); i++) {
        // Convert current and rest positions to SE3 transformations
        liegroups::SE3d T_current, T_rest;
        
        // Current configuration
        T_current.translation() = Vec3d(x[i][0], x[i][1], x[i][2]);
        liegroups::SO3d R_current = liegroups::SO3d::exp(Vec3d(x[i][3], x[i][4], x[i][5]));
        T_current.rotation() = R_current.matrix();
        
        // Rest configuration
        T_rest.translation() = Vec3d(x0[i][0], x0[i][1], x0[i][2]);
        liegroups::SO3d R_rest = liegroups::SO3d::exp(Vec3d(x0[i][3], x0[i][4], x0[i][5]));
        T_rest.rotation() = R_rest.matrix();
        
        // Compute relative transformation
        liegroups::SE3d T_rel = T_current * T_rest.inverse();
        
        // Get twist coordinates and compute force
        Deriv twist = T_rel.log();
        f[i] = -(this->m_K_section66 * twist) * lengths[i];
    }

    d_f.endEdit();
}

template<>
void BeamHookeLawForceField<defaulttype::Vec6Types>::addDForce(
    const MechanicalParams* mparams,
    DataVecDeriv& d_df,
    const DataVecDeriv& d_dx)
{
    WriteAccessor<DataVecDeriv> df = d_df;
    ReadAccessor<DataVecDeriv> dx = d_dx;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    const vector<Real>& lengths = this->d_length.getValue();

    df.resize(dx.size());
    
    // Linear approximation for small displacements
    for (size_t i = 0; i < dx.size(); i++) {
        df[i] = -(this->m_K_section66 * dx[i]) * kFactor * lengths[i];
    }
}

template<>
void BeamHookeLawForceField<defaulttype::Vec6Types>::addKToMatrix(
    const MechanicalParams* mparams,
    const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    const vector<Real>& lengths = this->d_length.getValue();

    const VecCoord& pos = this->mstate->read(core::vec_id::read_access::position)->getValue();
    
    // Add stiffness contribution for each beam element
    for (size_t n = 0; n < pos.size(); n++) {
        for (unsigned int i = 0; i < 6; i++) {
            for (unsigned int j = 0; j < 6; j++) {
                mat->add(offset + i + 6 * n, offset + j + 6 * n,
                        -kFact * this->m_K_section66[i][j] * lengths[n]);
            }
        }
    }
}

template<>
double BeamHookeLawForceField<defaulttype::Vec6Types>::getPotentialEnergy(
    const MechanicalParams* mparams,
    const DataVecCoord& x) const
{
    SOFA_UNUSED(mparams);
    
    const VecCoord& pos = x.getValue();
    const VecCoord& rest_pos = this->mstate->read(core::vec_id::read_access::restPosition)->getValue();
    const vector<Real>& lengths = this->d_length.getValue();
    
    double energy = 0.0;
    
    // Compute potential energy using LieGroups
    for (size_t i = 0; i < pos.size(); i++) {
        // Convert configurations to SE3
        liegroups::SE3d T_current, T_rest;
        
        // Current configuration
        T_current.translation() = Vec3d(pos[i][0], pos[i][1], pos[i][2]);
        liegroups::SO3d R_current = liegroups::SO3d::exp(Vec3d(pos[i][3], pos[i][4], pos[i][5]));
        T_current.rotation() = R_current.matrix();
        
        // Rest configuration
        T_rest.translation() = Vec3d(rest_pos[i][0], rest_pos[i][1], rest_pos[i][2]);
        liegroups::SO3d R_rest = liegroups::SO3d::exp(Vec3d(rest_pos[i][3], rest_pos[i][4], rest_pos[i][5]));
        T_rest.rotation() = R_rest.matrix();
        
        // Compute relative transformation
        liegroups::SE3d T_rel = T_current * T_rest.inverse();
        
        // Get twist coordinates
        Deriv twist = T_rel.log();
        
        // Compute energy contribution
        energy += 0.5 * (twist * (this->m_K_section66 * twist)) * lengths[i];
    }
    
    return energy;
}

using namespace sofa::defaulttype;

template class BeamHookeLawForceField<Vec3Types>;
template class BeamHookeLawForceField<Vec6Types>;

} // namespace sofa::component::forcefield

namespace Cosserat
{

void registerBeamHookeLawForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(
        sofa::core::ObjectRegistrationData(
            "This component is used to compute internal stress for torsion (along x) and bending (y and z)")
            .add<sofa::component::forcefield::BeamHookeLawForceField<sofa::defaulttype::Vec3Types>>(true)
            .add<sofa::component::forcefield::BeamHookeLawForceField<sofa::defaulttype::Vec6Types>>());
}

} // namespace Cosserat
