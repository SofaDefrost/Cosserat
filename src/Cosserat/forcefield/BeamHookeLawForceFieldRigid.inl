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
#include "BeamHookeLawForceFieldRigid.h"
#include <Cosserat/forcefield/BeamHookeLawForceFieldRigid.h>

#include <sofa/linearalgebra/FullVector.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/OptionsGroup.h> // ??
#include <iostream>
#include <algorithm>
#include <ctime>

using sofa::core::behavior::MechanicalState ;
using sofa::core::objectmodel::BaseContext ;
using sofa::helper::ReadAccessor ;
using sofa::helper::WriteAccessor ;
using sofa::core::VecCoordId;
using std::cout ;
using std::endl ;

namespace sofa::component::forcefield
{
    using sofa::core::behavior::MultiMatrixAccessor ;
    using sofa::core::behavior::BaseMechanicalState ;
    using sofa::helper::WriteAccessor ;


    template<typename DataTypes>
    BeamHookeLawForceFieldRigid<DataTypes>::BeamHookeLawForceFieldRigid()
            : Inherit1(),
              d_crossSectionShape( initData(&d_crossSectionShape, {"circular","rectangular"},
                                            "crossSectionShape",
                                            "shape of the cross-section. Can be: circular (tube with external radius being radius and internal radius being innerRadius ) or rectangular (lengthY and lengthZ) . Default is circular" )),
              d_youngModulus( initData( &d_youngModulus, 1.0e9, "youngModulus", "Young Modulus describes the stiffness of the material")),
              d_poissonRatio( initData( &d_poissonRatio, 0.45, "poissonRatio", "poisson Ratio describes the compressibility of the material")),
              d_length( initData( &d_length, "length", "length of each beam")),
              d_radius( initData( &d_radius, 1.0, "radius", "external radius of the cross section (if circular)")),
              d_innerRadius( initData( &d_innerRadius, 0.0, "innerRadius", "internal radius of the cross section (if circular)")),
              d_lengthY( initData( &d_lengthY, 1.0, "lengthY", "side length of the cross section along local y axis (if rectangular)")),
              d_lengthZ( initData( &d_lengthZ, 1.0, "lengthZ", "side length of the cross section along local z axis (if rectangular)")),
              d_variantSections(initData(&d_variantSections, false, "variantSections", "In case we have variant beam sections this has to be set to true")),
              d_youngModulusList(initData(&d_youngModulusList, "youngModulusList", "The list of Young modulus in case we have sections with variable physical properties")),
              d_poissonRatioList(initData(&d_poissonRatioList, "poissonRatioList", "The list of poisson's ratio in case we have sections with variable physical properties")),
              d_useInertiaParams(initData(&d_useInertiaParams, false, "useInertiaParams", "If the inertia parameters are given by the user, there is no longer any need to use @d_youngModulus and @d_poissonRatio.")),
              d_GI(initData(&d_GI, "GI", "The inertia parameter, GI")),
              d_GA(initData(&d_GA, "GA", "The inertia parameter, GA")),
              d_EA(initData(&d_EA, "EA", "The inertia parameter, EA")),
              d_EIy(initData(&d_EIy, "EIy", "The inertia parameter, EIy")),
               d_EIz(initData(&d_EIz, "EIz", "The inertia parameter, EIz")),
              d_buildTorsion(initData(&d_buildTorsion, true,"build_torsion", "build torsion or the elongation of the beam ?"))

    {
        compute_df=true;
    }

    template<typename DataTypes>
    BeamHookeLawForceFieldRigid<DataTypes>::~BeamHookeLawForceFieldRigid()
    {}

    template<typename DataTypes>
    void BeamHookeLawForceFieldRigid<DataTypes>::init()
    {
        Inherit1::init();
        reinit();
    }

    /*Cross-Section Properties Initialization: The reinit function begins by recalculating the properties
    related to the cross-section of the beams. It calculates the area moment of inertia (Iy and Iz),
    the polar moment of inertia (J), and the cross-sectional area (A).
    These calculations depend on the chosen cross-section shape, either circular or rectangular.
    The formulas used for these calculations are based on standard equations for these properties.*/
    template<typename DataTypes>
    void BeamHookeLawForceFieldRigid<DataTypes>::reinit()
    {
        // Precompute and store values
        Real Iy, Iz, J, A;
        if ( d_crossSectionShape.getValue().getSelectedItem() == "rectangular")  //rectangular cross-section
        {
            Real Ly = d_lengthY.getValue();
            Real Lz = d_lengthZ.getValue();

            const Real LyLzLzLz = Ly * Lz * Lz * Lz;
            const Real LzLyLyLy = Lz * Ly * Ly * Ly;

            Iy = LyLzLzLz / 12.0;
            Iz = LzLyLyLy / 12.0;
            J  = Iy + Iz;
            A  = Ly * Lz;

        }
        else //circular cross-section
        {
            msg_info() << "Cross section shape." << d_crossSectionShape.getValue().getSelectedItem() ;

            Real r = d_radius.getValue();
            Real rInner = d_innerRadius.getValue();
            const Real r4 = r * r * r * r;
            const Real rInner4 = rInner * rInner * rInner * rInner;

            Iy = M_PI * (r4 - rInner4) / 4.0;
            Iz = Iy;
            J = Iy + Iz;
            A = M_PI * (r * r - rInner4);

        }
        m_crossSectionArea = A;

        if( d_useInertiaParams.getValue() )
        {
            msg_info("BeamHookeLawForceFieldRigid")<< "Pre-calculated inertia parameters are used for the computation of the stiffness matrix.";
            m_K_section66[0][0] = d_GI.getValue();
            m_K_section66[1][1] = d_EIy.getValue();
            m_K_section66[2][2] = d_EIz.getValue();
            m_K_section66[3][3] = d_EA.getValue();
            m_K_section66[4][4] = d_GA.getValue();
            m_K_section66[5][5] = d_GA.getValue();
        }
        else
        {
            //Pour du vec 6, on a  _m_K =diag([G*J E*Iy E*Iz E*A G*As G*As]); % stifness matrix
            Real E = d_youngModulus.getValue();
            Real G = E/(2.0*(1.0+d_poissonRatio.getValue()));
            //Translational Stiffness (X, Y, Z directions):
            // Gs * J(i): Shearing modulus times the second moment of the area (torsional stiffness).
            // E * Js(i): Young's modulus times the second moment of the area (bending stiffness).
            m_K_section66[0][0] = G*J;
            m_K_section66[1][1] = E*Iy;
            m_K_section66[2][2] = E*Iz;
            //Rotational Stiffness (X, Y, Z directions):
            //E * A: Young's modulus times the cross-sectional area (axial stiffness).
            //Gs * A: Shearing modulus times the cross-sectional area (shearing stiffness).
            m_K_section66[3][3] = E*A;
            m_K_section66[4][4] = G*A;
            m_K_section66[5][5] = G*A;
        }
    }

    template<typename DataTypes>
    void BeamHookeLawForceFieldRigid<DataTypes>::addForce(const MechanicalParams* mparams,
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
        unsigned int sz = d_length.getValue().size();
        if(x.size()!= sz){
            msg_warning("BeamHookeLawForceField")<<" length : "<< sz <<"should have the same size as x... "<< x.size() <<"\n";
            compute_df = false;
            return;
        }
        for (unsigned int i=0; i<x.size(); i++)
            f[i] -= (m_K_section66 * (x[i] - x0[i])) * d_length.getValue()[i];

        d_f.endEdit();

    }

    template<typename DataTypes>
    void BeamHookeLawForceFieldRigid<DataTypes>::addDForce(const MechanicalParams* mparams,
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
            df[i] -= (m_K_section66 * dx[i])*kFactor* d_length.getValue()[i];
    }

    template<typename DataTypes>
    void BeamHookeLawForceFieldRigid<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                              const MultiMatrixAccessor* matrix)
    {
        MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
        BaseMatrix* mat = mref.matrix;
        unsigned int offset = mref.offset;
        Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

        const VecCoord& pos = this->mstate->read(core::ConstVecCoordId::position())->getValue();
        for (unsigned int n=0; n<pos.size(); n++)
        {
            for(unsigned int i = 0; i < 6; i++)
                for (unsigned int j = 0; j< 6; j++)
                    mat->add(offset + i + 6*n, offset + j + 6*n, -kFact * m_K_section66[i][j]*d_length.getValue()[n]);

        }
    }

    template<typename DataTypes>
    double BeamHookeLawForceFieldRigid<DataTypes>::getPotentialEnergy(const MechanicalParams* mparams,
                                                                 const DataVecCoord& d_x) const
    {
        SOFA_UNUSED(mparams);
        SOFA_UNUSED(d_x);

        return 0.0;
    }
    template<typename DataTypes>
    typename BeamHookeLawForceFieldRigid<DataTypes>::Real BeamHookeLawForceFieldRigid<DataTypes>::getRadius()
    {
        return d_radius.getValue();
    }

} // forcefield
