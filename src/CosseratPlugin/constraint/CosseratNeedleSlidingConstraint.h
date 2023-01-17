/******************************************************************************
 *               SOFA, Simulation Open-Framework Architecture                  *
 *                (c) 2020 INRIA, USTL, UJF, CNRS, MGH                    *
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
 *                      Plugin Cosserat v1.0                                   *
 *                                                                             *
 * This plugin is also distributed under the GNU LGPL (Lesser General          *
 * Public License) license with the same conditions than SOFA.                 *
 *                                                                             *
 * Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
 *               Ecole Centrale de Lille)                                      *
 *                                                                             *
 * Contact information: https://project.inria.fr/softrobot/contact/            *
 *                     adagolodjo@protonmail.com                               *
 ******************************************************************************/

#pragma once

#include "../config.h.in"

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/helper/OptionsGroup.h>
#include <sofa/core/behavior/Constraint.h>

namespace sofa::component::constraintset
{

using sofa::core::ConstraintParams;
using sofa::core::VecCoordId;
using sofa::core::behavior::Constraint;
using sofa::core::objectmodel::Data;
using sofa::core::visual::VisualParams;
using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec3fTypes;
using sofa::helper::ReadAccessor;
using sofa::linearalgebra::BaseVector;
using sofa::core::MultiVecDerivId ;

using sofa::defaulttype::Rigid3Types ;
using sofa::defaulttype::Vec3Types ;


template<class T>
class CosseratNeedleSlidingConstraintSpecialization {};


/**
 * This class contains common implementation of cable constraints
*/
template< class DataTypes >
class CosseratNeedleSlidingConstraint : public Constraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CosseratNeedleSlidingConstraint,DataTypes), SOFA_TEMPLATE(Constraint,DataTypes));

    /// That any templates variation of CosseratNeedleSlidingConstraintSpecialization are friend.
    template<typename>
    friend class CosseratNeedleSlidingConstraintSpecialization ;

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename Coord::value_type Real;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;

    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    typedef Data<MatrixDeriv> DataMatrixDeriv;
    typedef type::vector<unsigned int> SetIndexArray;


public:
    CosseratNeedleSlidingConstraint(MechanicalState* object = nullptr);

    ~CosseratNeedleSlidingConstraint() override;

    ////////////////////////// Inherited from BaseObject ////////////////////
    void init() override;
    void reinit() override;
    void draw(const VisualParams* vparams) override;
    /////////////////////////////////////////////////////////////////////////

    ////////////////////////// Inherited from Actuator //////////////////////
    void buildConstraintMatrix(const ConstraintParams* cParams,
                               DataMatrixDeriv &cMatrix,
                               unsigned int &cIndex,
                               const DataVecCoord &x) override;
    void getConstraintViolation(const ConstraintParams* cParams,
                                BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &v) override;
//    void getConstraintViolation(const ConstraintParams* cParams,
//                                BaseVector *resV, const DataVecCoord &x, const DataVecDeriv &v) override;
    void getConstraintResolution(const ConstraintParams*,
                                 std::vector<core::behavior::ConstraintResolution*>& resTab,
                                 unsigned int& offset) override;

    /////////////////////////////////////////////////////////////////////////

    ////////////////////////// Inherited from BaseConstraint ////////////////
    //    void storeLambda(const ConstraintParams* cParams,
    //                     core::MultiVecDerivId res,
    //                     const BaseVector* lambda) override;
    /////////////////////////////////////////////////////////////////////////

protected:
    //Input data
    Data<type::vector< Real > >   d_value;
    Data<unsigned int>            d_valueIndex;
    Data<helper::OptionsGroup>    d_valueType;
    // This data fields indicates which DoFs are taken into account when applying
    // the constraint. In case of nonRigid template instanciation (typically Vec3),
    // only the first 3 boolean are used.
    Data<type::Vec<6,bool>>       d_useDirections;
    // displacement = the constraint will impose the displacement provided in data d_inputValue[d_iputIndex]
    // force = the constraint will impose the force provided in data d_inputValue[d_iputIndex]


protected:

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using Constraint<DataTypes>::d_componentState ;
    ////////////////////////////////////////////////////////////////////////////
    /// \brief internalInit
    unsigned int m_nbLines ;
    unsigned int m_constraintId ;

    void internalInit();
};

// Declares template as extern to avoid the code generation of the template for
// each compilation unit. see: http://www.stroustrup.com/C++11FAQ.html#extern-templates
//extern template class CosseratNeedleSlidingConstraint<defaulttype::Vec3Types>;
#if !defined(SOFA_COMPONENT_CONSTRAINTSET_COSSERATNEEDLESLIDINGCONSTRAINT_CPP)
extern template class SOFA_COSSERATPLUGIN_API CosseratNeedleSlidingConstraint< Vec3Types >;
extern template class SOFA_COSSERATPLUGIN_API CosseratNeedleSlidingConstraint< Rigid3Types >;
#endif

} // namespace sofa
