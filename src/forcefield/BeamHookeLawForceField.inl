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
*                           Plugin SoftRobots    v1.0                         *
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

#include "BeamHookeLawForceField.h"
#include <SofaBaseLinearSolver/FullVector.h>
#include <sofa/core/behavior/MechanicalState.h>

using sofa::core::behavior::MechanicalState ;
using sofa::core::objectmodel::BaseContext ;

#include <sofa/core/behavior/ForceField.inl>
using sofa::helper::ReadAccessor ;
using sofa::helper::WriteAccessor ;
using sofa::core::VecCoordId;

// ??
#include <sofa/helper/OptionsGroup.h>


#include <iostream>
using std::cout ;
using std::endl ;

#include <algorithm>
#include <ctime>

namespace sofa::component::forcefield
{

using sofa::component::linearsolver::DefaultMultiMatrixAccessor ;
using sofa::core::behavior::MultiMatrixAccessor ;
using sofa::core::behavior::BaseMechanicalState ;
using sofa::helper::WriteAccessor ;


template<typename DataTypes>
BeamHookeLawForceField<DataTypes>::BeamHookeLawForceField()
    : Inherit1(),
      d_crossSectionShape( initData(&d_crossSectionShape, OptionsGroup(2,"circular","rectangular"),
                                    "crossSectionShape",
                                    "shape of the cross-section. Can be: circular (tube with external radius being radius and internal radius being innerRadius ) or rectangular (lengthY and lengthZ) . Default is circular" )),
      d_youngModulus( initData( &d_youngModulus, 1.0e9, "youngModulus", "Young Modulus describes the stiffness of the material")),
      d_poissonRatio( initData( &d_poissonRatio, 0.45, "poissonRatio", "poisson Ratio describes the compressibility of the material")),
      d_length( initData( &d_length, "length", "lenght of each beam")),
      d_radius( initData( &d_radius, 1.0, "radius", "external radius of the cross section (if circular)")),
      d_innerRadius( initData( &d_innerRadius, 0.0, "innerRadius", "internal radius of the cross section (if circular)")),
      d_lengthY( initData( &d_lengthY, 1.0, "lengthY", "side length of the cross section along local y axis (if rectangular)")),
      d_lengthZ( initData( &d_lengthZ, 1.0, "lengthZ", "side length of the cross section along local z axis (if rectangular)")),
      d_varianteSections( initData( &d_varianteSections, false, "varianteSections", "In case we have variante beam section this has to be set to true")),
      d_youngModulusList(initData(&d_youngModulusList, "youngModulusList", "The list of Young modulus in case we have sections with variable physical properties")),
      d_poissonRatioList(initData(&d_poissonRatioList, "poissonRatioList", "The list of poisson's ratio in case we have sections with variable physical properties"))
{
    compute_df=true;
}


template<typename DataTypes>
BeamHookeLawForceField<DataTypes>::~BeamHookeLawForceField()
{}

template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::init()
{
    Inherit1::init();

    reinit();
}

template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::reinit()
{
    Real Iy, Iz, J;
    Real A;
    if ( d_crossSectionShape.getValue().getSelectedItem() == "rectangular" )
    {
        Real Ly = d_lengthY.getValue();
        Real Lz = d_lengthZ.getValue();

        Iy=Ly*Lz*Lz*Lz/12.0;
        Iz=Lz*Ly*Ly*Ly/12.0;
        J=Iy + Iz;
        A = Ly*Lz;

    }
    else //circular section
    {
        msg_info() << "Cross section shape." << d_crossSectionShape.getValue().getSelectedItem() ;

        Real r = d_radius.getValue();
        Real rInner = d_innerRadius.getValue();
        Iz = M_PI*(r*r*r*r - rInner*rInner*rInner*rInner)/4.0;
        Iy = Iz ;
        J = Iz + Iy;
        A = M_PI*(r*r - rInner*rInner);

    }
    m_crossSectionArea = A;

    if(!d_varianteSections.getValue()){

        Real E= d_youngModulus.getValue();
        Real G= E/(2.0*(1.0+d_poissonRatio.getValue()));

        m_K_section[0][0] = G*J;
        m_K_section[1][1] = E*Iy;
        m_K_section[2][2] = E*Iz;
    }else {
        msg_info("BeamHookeLawForceField")<< "=====> Multi section";
        m_K_sectionList.clear();
        size_t szYM = d_youngModulusList.getValue().size();
        size_t szPR = d_poissonRatioList.getValue().size();
        size_t szL  = d_length.getValue().size();

        std::cout<< "BeamHookeLawForceField :" << "szYM : "<< szYM << " szPR "<< szPR << "  szL " << szL<< std::endl;

        if((szL != szPR)||(szL != szYM)){
            msg_error("BeamHookeLawForceField")<< "Please, lenght, youngModulusList and poissonRatioList should have the same size";
            return;
        }

        for(size_t k=0; k<szL; k++){
            Mat33 _m_K_section;
            Real E= d_youngModulusList.getValue()[k];
            Real G= E/(2.0*(1.0+d_poissonRatioList.getValue()[k]));

            _m_K_section[0][0] = G*J;
            _m_K_section[1][1] = E*Iy;
            _m_K_section[2][2] = E*Iz;
            m_K_sectionList.push_back(_m_K_section);
        }
    }
}

template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::addForce(const MechanicalParams* mparams,
                                                 DataVecDeriv& d_f,
                                                 const DataVecCoord& d_x,
                                                 const DataVecDeriv& d_v)
{
    SOFA_UNUSED(d_v);
    SOFA_UNUSED(mparams);

    if(!this->getMState()) {
        msg_info("BeamHookeLawForceField") << "No Mechanical State found, no force will be computed..." << "\n";
        compute_df=false;
        return;
    }
    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    // get the rest position (for non straight shape)
    const VecCoord& x0 = this->mstate->read(VecCoordId::restPosition())->getValue();

    f.resize(x.size());
    if(x.size()!=d_length.getValue().size()){
        msg_warning("BeamHookeLawForceField")<<" length should have the same size as x..."<<"\n";
        compute_df = false;
        return;
    }

    if(!d_varianteSections.getValue())
        for (unsigned int i=0; i<x.size(); i++)
            f[i] -= (m_K_section * (x[i] - x0[i])) * d_length.getValue()[i];
    else
        for (unsigned int i=0; i<x.size(); i++)
            f[i] -= (m_K_sectionList[i] * (x[i] - x0[i])) * d_length.getValue()[i];
    d_f.endEdit();

}

template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::addDForce(const MechanicalParams* mparams,
                                                  DataVecDeriv&  d_df ,
                                                  const DataVecDeriv&  d_dx)
{
    if (!compute_df)
        return;

    WriteAccessor< DataVecDeriv > df = d_df;
    ReadAccessor< DataVecDeriv > dx = d_dx;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    df.resize(dx.size());
    if(!d_varianteSections.getValue())
        for (unsigned int i=0; i<dx.size(); i++)
            df[i] -= (m_K_section * dx[i])*kFactor* d_length.getValue()[i];
    else
        for (unsigned int i=0; i<dx.size(); i++)
            df[i] -= (m_K_sectionList[i] * dx[i])*kFactor* d_length.getValue()[i];
}

template<typename DataTypes>
double BeamHookeLawForceField<DataTypes>::getPotentialEnergy(const MechanicalParams* mparams,
                                                             const DataVecCoord& d_x) const
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_x);

    return 0.0;
}

template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                     const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const VecCoord& pos = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    for (unsigned int n=0; n<pos.size(); n++)
    {
        if(!d_varianteSections.getValue())
            for(int i = 0; i < 3; i++)
                for (int j = 0; j< 3; j++)
                    mat->add(offset + i + 3*n, offset + j + 3*n, -kFact * m_K_section[i][j]*d_length.getValue()[n]);
        else
            for(int i = 0; i < 3; i++)
                for (int j = 0; j< 3; j++)
                    mat->add(offset + i + 3*n, offset + j + 3*n, -kFact * m_K_sectionList[n][i][j] * d_length.getValue()[n]);
    }
}




} // forcefield
