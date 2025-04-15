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
#include <Cosserat/config.h>
#include <Cosserat/forcefield/BaseBeamHookeLawForceField.h>

#include <sofa/linearalgebra/CompressedRowSparseMatrix.h>
#include <sofa/helper/ScopedAdvancedTimer.h>

namespace sofa::component::forcefield
{

using sofa::type::Vec;
using sofa::type::Mat;
using sofa::type::vector;
using sofa::core::MechanicalParams;
using sofa::linearalgebra::BaseMatrix;
using sofa::linearalgebra::CompressedRowSparseMatrix;
using sofa::core::behavior::MultiMatrixAccessor;
/**
 * This component is used to compute the Hooke's law on a beam computed on strain / stress
 * Only bending and torsion strain / stress are considered here
 * It derives from BaseBeamHookeLawForceField to utilize Lie Group operations
*/
template<typename DataTypes>
class BeamHookeLawForceField : public BaseBeamHookeLawForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BeamHookeLawForceField, DataTypes), SOFA_TEMPLATE(BaseBeamHookeLawForceField, DataTypes));
    
    typedef BaseBeamHookeLawForceField<DataTypes> Inherit1;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    
    typedef typename Inherit1::Vector Vector;
    typedef typename Inherit1::Vector3 Vector3;
    typedef typename Inherit1::SO3Type SO3Type;
    
    typedef CompressedRowSparseMatrix<Mat<3, 3, Real>> CSRMat33B66;
    
    typedef typename CompressedRowSparseMatrix<Mat<3, 3, Real>>::ColBlockConstIterator _3_3_ColBlockConstIterator;
    typedef typename CompressedRowSparseMatrix<Mat<3, 3, Real>>::RowBlockConstIterator _3_3_RowBlockConstIterator;
    typedef typename CompressedRowSparseMatrix<Mat<3, 3, Real>>::BlockConstAccessor _3_3_BlockConstAccessor;
    typedef typename CompressedRowSparseMatrix<Mat<3, 3, Real>>::BlockAccessor _3_3_BlockAccessor;

protected:
    // Implementation of abstract methods from base class
    Vector getPosition(const Coord& coord) const override;
    SO3Type getRotation(const Coord& coord) const override;
    Vector getForce(const Deriv& deriv) const override;
    Vector getMoment(const Deriv& deriv) const override;
    Deriv createDeriv(const Vector& force, const Vector& moment) const override;

public:
    BeamHookeLawForceField();
    virtual ~BeamHookeLawForceField() = default;
private:
    // No additional member variables or methods as we use the ones from BaseBeamHookeLawForceField

private :

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using ForceField<DataTypes>::getContext ;
    using ForceField<DataTypes>::f_printLog ;
    ////////////////////////////////////////////////////////////////////////////
};

#if !defined(SOFA_COSSERAT_CPP_BeamHookeLawForceField)
  extern template class SOFA_COSSERAT_API BeamHookeLawForceField<defaulttype::Vec3Types>;
  extern template class SOFA_COSSERAT_API BeamHookeLawForceField<defaulttype::Vec6Types>;
#endif

} // forcefield
