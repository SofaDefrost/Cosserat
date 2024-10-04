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

#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/linearalgebra/BaseMatrix.h>
#include <sofa/linearalgebra/CompressedRowSparseMatrix.h>
#include <sofa/type/Mat.h>
#include <sofa/type/Vec.h>

namespace sofa::component::forcefield {

using sofa::core::MechanicalParams;
using sofa::core::behavior::ForceField;
using sofa::core::behavior::MultiMatrixAccessor;
using sofa::linearalgebra::BaseMatrix;
using sofa::linearalgebra::CompressedRowSparseMatrix;
using sofa::type::Mat;
using sofa::type::Vec;
using sofa::type::vector;

using sofa::helper::OptionsGroup;

/**
 * This component is used to compute the Hooke's law on a beam computed on
 * strain / stress Only bending and torsion strain / stress are considered here
 */
template <typename DataTypes>
class BeamHookeLawForceField : public ForceField<DataTypes> {
public:
  SOFA_CLASS(SOFA_TEMPLATE(BeamHookeLawForceField, DataTypes),
             SOFA_TEMPLATE(ForceField, DataTypes));

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

  typedef CompressedRowSparseMatrix<Mat33> CSRMat33B66;

  typedef typename CompressedRowSparseMatrix<Mat33>::ColBlockConstIterator
      _3_3_ColBlockConstIterator;
  typedef typename CompressedRowSparseMatrix<Mat33>::RowBlockConstIterator
      _3_3_RowBlockConstIterator;
  typedef typename CompressedRowSparseMatrix<Mat33>::BlockConstAccessor
      _3_3_BlockConstAccessor;
  typedef typename CompressedRowSparseMatrix<Mat33>::BlockAccessor
      _3_3_BlockAccessor;

public:
  BeamHookeLawForceField();
  virtual ~BeamHookeLawForceField();

  ////////////////////////// Inherited from BaseObject /////////////////////////
  void init() override;
  void reinit() override;
  ///////////////////////////////////////////////////////////////////////////

  ////////////////////////// Inherited from ForceField /////////////////////////
  void addForce(const MechanicalParams *mparams, DataVecDeriv &f,
                const DataVecCoord &x, const DataVecDeriv &v) override;

  void addDForce(const MechanicalParams *mparams, DataVecDeriv &df,
                 const DataVecDeriv &dx) override;

  void addKToMatrix(const MechanicalParams *mparams,
                    const MultiMatrixAccessor *matrix) override;

  double getPotentialEnergy(const MechanicalParams *mparams,
                            const DataVecCoord &x) const override;
  ////////////////////////////////////////////////////////////////////////////

  Real getRadius();

protected:
  Data<helper::OptionsGroup> d_crossSectionShape;
  Data<Real> d_youngModulus;         /// youngModulus
  Data<Real> d_poissonRatio;         /// poissonRatio
  Data<type::vector<Real>> d_length; /// length of each beam

  /// Circular Cross Section
  Data<Real> d_radius;
  Data<Real> d_innerRadius;
  /// Rectangular Cross Section
  Data<Real> d_lengthY;
  Data<Real> d_lengthZ;
  // In case we have a beam with different properties per section
  Data<bool> d_variantSections; /// bool to identify different beams sections
  Data<type::vector<Real>> d_youngModulusList; /// youngModulus
  Data<type::vector<Real>> d_poissonRatioList; /// poissonRatio
  /// If the inertia parameters are given by the user, there is no longer any
  /// need to use YG.
  Data<bool> d_useInertiaParams;
  Data<Real> d_GI;
  Data<Real> d_GA;
  Data<Real> d_EA;
  Data<Real> d_EI;
  Data<Real> d_EIy;
  Data<Real> d_EIz;
  Data<Real> d_MomQuad;

  bool compute_df;
  Mat33 m_K_section;
  Mat66 m_K_section66;
  type::vector<Mat33> m_K_sectionList;

  /// Cross-section area
  Real m_crossSectionArea;

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

#if !defined(SOFA_COSSERAT_CPP_BeamHookeLawForceField)
extern template class SOFA_COSSERAT_API
    BeamHookeLawForceField<defaulttype::Vec3Types>;
extern template class SOFA_COSSERAT_API
    BeamHookeLawForceField<defaulttype::Vec6Types>;
#endif

} // namespace sofa::component::forcefield
