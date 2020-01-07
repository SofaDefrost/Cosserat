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
#ifndef SOFA_COMPONENT_FORCEFIELD_CosseratInternalActuation_H
#define SOFA_COMPONENT_FORCEFIELD_CosseratInternalActuation_H

#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/MechanicalParams.h>
#include <SofaBaseLinearSolver/CompressedRowSparseMatrix.h>
#include <SofaBaseLinearSolver/DefaultMultiMatrixAccessor.h>

#include <sofa/helper/OptionsGroup.h>

namespace sofa
{

namespace component
{

namespace forcefield
{
using sofa::defaulttype::Vec ;
using sofa::defaulttype::Mat ;
using sofa::helper::vector;
using sofa::core::MechanicalParams;
using sofa::defaulttype::BaseMatrix;
using sofa::core::behavior::ForceField ;
using sofa::component::linearsolver::CompressedRowSparseMatrix ;
using sofa::core::behavior::MultiMatrixAccessor ;

using sofa::helper::OptionsGroup;

/**
 * This component is used to compute the Hooke's law on a beam computed on strain / stress
 * Only bending and torsion strain / stress are considered here
*/
template<typename DataTypes>
class CosseratInternalActuation : public ForceField<DataTypes>
{
public :
    SOFA_CLASS(SOFA_TEMPLATE(CosseratInternalActuation, DataTypes), SOFA_TEMPLATE(ForceField, DataTypes));

    typedef typename DataTypes::Real     Real;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord;
    typedef typename DataTypes::Deriv    Deriv;

    typedef Data<VecCoord>    DataVecCoord;
    typedef Data<VecDeriv>    DataVecDeriv;

    typedef Vec<3, Real>                Vec3;
    typedef Mat<3, 3, Real>             Mat33;

    typedef CompressedRowSparseMatrix<Mat33> CSRMat33B66;
    typedef defaulttype::Vector3 vector3;

    typedef typename CompressedRowSparseMatrix<Mat33>::ColBlockConstIterator _3_3_ColBlockConstIterator;
    typedef typename CompressedRowSparseMatrix<Mat33>::RowBlockConstIterator _3_3_RowBlockConstIterator;
    typedef typename CompressedRowSparseMatrix<Mat33>::BlockConstAccessor _3_3_BlockConstAccessor;
    typedef typename CompressedRowSparseMatrix<Mat33>::BlockAccessor _3_3_BlockAccessor;


public :
    CosseratInternalActuation();
    virtual ~CosseratInternalActuation();

    ////////////////////////// Inherited from BaseObject /////////////////////////
    void init() override;
    ///////////////////////////////////////////////////////////////////////////

    //    void computeArgument(const VecCoord &distance, const VecCoord &x, const int id, Vec3 &argu, const double &C);
    void computeIntegrale(const double &Li, const double &Li_1, const VecCoord &x, const int id, Coord & integral);

    ////////////////////////// Inherited from ForceField /////////////////////////
    void addForce(const MechanicalParams* mparams,
                  DataVecDeriv& f ,
                  const DataVecCoord& x ,
                  const DataVecDeriv& v) override;

    void addDForce(const MechanicalParams* mparams,
                   DataVecDeriv&   df ,
                   const DataVecDeriv&
                   dx ) override;


    void addKToMatrix(const MechanicalParams* mparams,
                      const MultiMatrixAccessor* matrix) override;

    double getPotentialEnergy(const MechanicalParams* mparams,
                              const DataVecCoord& x) const override;
    ////////////////////////////////////////////////////////////////////////////

protected:
    Data<helper::OptionsGroup>   d_crossSectionShape;

    Data<double>                      d_youngModululs; /// youngModulus
    Data<double>                      d_poissonRatio; /// poissonRatio

    Data<helper::vector<double>>      d_length ; /// length of each beam

    /// Circular Cross Section
    Data<Real>          d_radius;
    Data<Real>          d_innerRadius;

    /// Rectangular Cross Section
    Data<Real>          d_lengthY;
    Data<Real>          d_lengthZ;

    //distance from the midlesection
    Data<helper::vector<Coord>>      d_distance0; // distance between the midleline and the calble
    Data<helper::vector<Coord>>      d_distance1; // distance between the midleline and the calble
    Data<helper::vector<Coord>>      d_ddistance0; // the derivative of the distance between the midleline and the calble with respect to x
    Data<helper::vector<Coord>>      d_ddistance1; // the derivative of the distance between the midleline and the calble with respect to x
    Data<double>                     d_Tt ; // Cable tension
    Data<helper::vector<Coord>>      d_integral; // the derivative of the distance between the midleline and the calble with respect to x


private :

    Mat33 m_K_section;
    bool compute_df;

    //Gaussian quadrature parameters for 2 points

    defaulttype::Vec2 m_gaussCoeff = defaulttype::Vec2(1.0/sqrt(3.0),0.57735); // Gauss quadrature coefficients
    defaulttype::Vec2 m_gaussWeights = defaulttype::Vec2(1.0,1.0);   //Gauss quadrature weights



    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using ForceField<DataTypes>::getContext ;
    using ForceField<DataTypes>::f_printLog ;
    using ForceField<DataTypes>::mstate ;
    ////////////////////////////////////////////////////////////////////////////
};


} // forcefield
} // component
} // sofa

#endif // SOFA_COMPONENT_FORCEFIELD_CosseratInternalActuation_H
