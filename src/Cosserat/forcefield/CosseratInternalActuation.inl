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
#include "CosseratInternalActuation.h"

#include <sofa/linearalgebra/FullVector.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/OptionsGroup.h>

#include <algorithm>
#include <ctime>
#include <iostream>

using sofa::core::behavior::MechanicalState ;
using sofa::core::objectmodel::BaseContext ;
using sofa::helper::ReadAccessor ;
using sofa::helper::WriteAccessor ;
using sofa::core::VecCoordId;
using std::cout ;
using std::endl ;

namespace sofa::component::forcefield
{

using sofa::core::behavior::DefaultMultiMatrixAccessor ;
using sofa::core::behavior::MultiMatrixAccessor ;
using sofa::core::behavior::BaseMechanicalState ;
using sofa::helper::WriteAccessor ;


template<typename DataTypes>
CosseratInternalActuation<DataTypes>::CosseratInternalActuation()
    : Inherit1(),
      d_crossSectionShape( initData(&d_crossSectionShape, {"circular","rectangular"},
                                    "crossSectionShape",
                                    "shape of the cross-section. Can be: circular (tube with external radius being radius and internal radius being innerRadius ) or rectangular (lengthY and lengthZ) . Default is circular" )),
      d_youngModululs( initData( &d_youngModululs, 1.0e6, "youngModulus", "Young Modulus describes the stiffness of the material")),
      d_poissonRatio( initData( &d_poissonRatio, 0.45, "poissonRatio", "poisson Ratio describes the compressibility of the material")),
      d_length( initData( &d_length, "length", "lenght of each beam")),
      d_radius( initData( &d_radius, 1.0, "radius", "external radius of the cross section (if circular)")),
      d_innerRadius( initData( &d_innerRadius, 0.0, "innerRadius", "internal radius of the cross section (if circular)")),
      d_lengthY( initData( &d_lengthY, 1.0, "lengthY", "side length of the cross section along local y axis (if rectangular)")),
      d_lengthZ( initData( &d_lengthZ, 1.0, "lengthZ", "side length of the cross section along local z axis (if rectangular)")),
      d_distance0( initData( &d_distance0,  "distance0", "distance between the midleline and the cable")),
      d_distance1( initData( &d_distance1,  "distance1", "distance between the midleline and the cable")),
      d_ddistance0( initData( &d_ddistance0,  "ddistance0", "the derivative of the distance between the midleline and the calble with respect to x")),
      d_ddistance1( initData( &d_ddistance1,  "ddistance1", "the derivative of the distance between the midleline and the calble with respect to x")),
      d_Tt( initData( &d_Tt,  "tension", "the cable tension according to t")),
      d_integral( initData( &d_integral,  "integral", "The value of the integral of all the tension"))
{
    compute_df=true;
}


template<typename DataTypes>
CosseratInternalActuation<DataTypes>::~CosseratInternalActuation()
{}

template<typename DataTypes>
void CosseratInternalActuation<DataTypes>::init()
{
    Inherit1::init();

    Real Iy, Iz, J, A;

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

    Real E= d_youngModululs.getValue();
    Real G= E/(2.0*(1.0+d_poissonRatio.getValue()));

    m_K_section[0][0] = G*J;
    m_K_section[1][1] = E*Iy;
    m_K_section[2][2] = E*Iz;


    //Compute the size of each beam
    //    helper::WriteAccessor<Data<helper::vector<double>>>  length = d_length;
    //    const VecCoord& pos = this->mstate->read(core::ConstVecCoordId::position())->getValue();
    //    helper::vector<double> m_vector; m_vector.clear();

    //    length.resize(pos.size());
    //    for (size_t i = 0; i < pos.size()-1; i++){
    //        if(i==0){
    //            length[i] = (pos[i]).norm();
    //            m_vector.push_back(0.0);
    //        }
    //        else{
    //            length[i] = (pos[i] - pos[i-1]).norm();
    //            m_vector.push_back(m_vector[i-1]+length[i]);
    //        }
    //    }

    //    std::cout << "============>length : \n"<< length << std::endl;
    //    std::cout << "============>distance : \n"<< m_vector << std::endl;

}

//template<typename DataTypes>
//void CosseratInternalActuation<DataTypes>::computeArgument(const VecCoord & distance, const VecCoord &x, const int id, Vec3 &argu, const double & C)
//{
//    //Compute s
//    //This computation is made manually and use for the computation of d(s) and d'(s)
//    //double s = ((Li-Li_1)/2.0)*C + (Li+Li_1)/2.0;

//    //    //compute argu
//    //    Coord ds = d_distance.getValue()[id];
//    //    Coord dds = d_ddistance.getValue()[id]; //derivative of the distance

//    //    //    std::cout << "x[id] :"<< x[id] << std::endl;
//    //    std::cout<< " The derivative is : "<<  dds << std::endl;
//    //    Coord vec = cross(x[id],ds) + Coord(1.0,0.0,0.0) + dds;
//    //    argu = cross(ds,vec)/(vec.norm());
//}

template<typename DataTypes>
void CosseratInternalActuation<DataTypes>::computeIntegrale(const double &Li, const double& Li_1, const VecCoord& x, const int id, Coord & integral)
{
    Coord arg0, arg1 ;
    VecCoord distance0 = d_distance0.getValue(); // d(X0) , X0 = [((Li-Li_1)/2)*C1 + (Li+Li_1)/2.0]
    VecCoord distance1 = d_distance1.getValue(); // d(X1) , X1 = [((Li-Li_1)/2)*C2 + (Li+Li_1)/2.0]

    VecCoord ddistance0 = d_ddistance0.getValue(); // the derivative of d
    VecCoord ddistance1 = d_ddistance1.getValue(); // the derivative of d

    //Compute f(X1) and f(X2) according to gauss parameters
    Coord vec0 = cross(x[id],distance0[id]) + Coord(1.0,0.0,0.0) + ddistance0[id];
    arg0 = cross(distance0[id],vec0)/(vec0.norm());
    //    std::cout << "fx0: "<< arg0 << std::endl;


    Coord vec1 = cross(x[id],distance1[id]) + Coord(1.0,0.0,0.0) + ddistance1[id];
    arg1 = cross(distance1[id],vec1)/(vec1.norm());
    //    std::cout << "fx1: "<< arg1 << std::endl;

    integral = ((Li-Li_1)/2.0) * (m_gaussWeights[0]*arg0 + m_gaussWeights[1]*arg1);

}

template<typename DataTypes>
void CosseratInternalActuation<DataTypes>::addForce(const MechanicalParams* mparams,
                                                    DataVecDeriv& d_f,
                                                    const DataVecCoord& d_x,
                                                    const DataVecDeriv& d_v)
{
    SOFA_UNUSED(d_v);
    SOFA_UNUSED(mparams);

    if(!this->getMState()) {
        msg_info("CosseratInternalActuation") << "No Mechanical State found, no force will be computed..." << "\n";
        compute_df=false;
        return;
    }
    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    // get the rest position (for non straight shape)
    const VecCoord& x0 = this->mstate->read(VecCoordId::restPosition())->getValue();

    f.resize(x.size());
    if(x.size()!=d_length.getValue().size()){
        msg_warning("CosseratInternalActuation")<<" length should have the same size as x..."<<"\n";
        compute_df = false;
        return;
    }

    for (unsigned int i=0; i<x.size(); i++)
    {
        //compute the tension internal force
        //(const double &Li, const double& Li_1, const VecCoord& x, const int id, Coord & integral)
        // Coord integral =  Coord(0.0,0.0,0.0);
        double Li = 0.0; double Li_1 = 0.0;

        for(unsigned j=0; j<=i; j++) Li += d_length.getValue()[j] ;
        if(i>0) for(unsigned j=0; j<i; j++) Li_1 += d_length.getValue()[j] ;

        //        computeIntegrale(Li,Li_1, x, i, integral);
        //        std::cout<< "Li_1 :"<< Li_1 << " ==> Li :"<< Li<<" ==> xi : "<< x[i]<< std::endl;
        //        std::cout << "Integral 0 :" << integral << std::endl;
        //        std::cout << "Integral 1 :" << d_integral.getValue()[i] << std::endl;
        //        f[i] -= (m_K_section * (x[i] - x0[i])) * d_length.getValue()[i] + d_Tt.getValue() * d_integral.getValue()[i];
        f[i] -= (m_K_section * (x[i] - x0[i])) * d_length.getValue()[i] ;
    }
    //    std::cout << "The finale force is : "<< f << std::endl;
    d_f.endEdit();

}

template<typename DataTypes>
void CosseratInternalActuation<DataTypes>::addDForce(const MechanicalParams* mparams,
                                                     DataVecDeriv&  d_df ,
                                                     const DataVecDeriv&  d_dx)
{
    if (!compute_df)
        return;

    WriteAccessor< DataVecDeriv > df = d_df;
    ReadAccessor< DataVecDeriv > dx = d_dx;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    df.resize(dx.size());
    for (unsigned int i=0; i<dx.size(); i++)
    {
        df[i] -= (m_K_section * dx[i])*kFactor* d_length.getValue()[i];
    }

}

template<typename DataTypes>
double CosseratInternalActuation<DataTypes>::getPotentialEnergy(const MechanicalParams* mparams,
                                                                const DataVecCoord& d_x) const
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_x);

    return 0.0;
}

template<typename DataTypes>
void CosseratInternalActuation<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                        const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    sofa::linearalgebra::BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const VecCoord& pos = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    for (unsigned int n=0; n<pos.size(); n++)
    {
        for(int i = 0; i < 3; i++)
        {
            for (int j = 0; j< 3; j++)
                mat->add(offset + i + 3*n, offset + j + 3*n, -kFact * m_K_section[i][j]*d_length.getValue()[n]);
        }
    }
}




} // forcefield
