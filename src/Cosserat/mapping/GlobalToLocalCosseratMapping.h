/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.        *
 *******************************************************************************
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
#pragma once
#include <Cosserat/config.h>
#include <Cosserat/mapping/BaseCosseratMapping.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cmath>
#include <sofa/core/BaseMapping.h>
#include <sofa/core/Mapping.h>
#include <sofa/core/State.h>
#include <sofa/helper/ColorMap.h>
#include <vector>

namespace Cosserat::mapping {
namespace {
using Mat3x6 = sofa::type::Mat<3, 6, SReal>;
using Mat6x3 = sofa::type::Mat<6, 3, SReal>;
using sofa::Data;
using sofa::type::Mat4x4;
using sofa::type::Mat6x6;
using sofa::type::Quat;
using sofa::type::Vec3;
using sofa::type::Vec6;
} // namespace

template <class TIn, class TOut>
class GlobalToLocalCosseratMapping : public sofa::core::Mapping<TIn, TOut> {
public:
  SOFA_CLASS(SOFA_TEMPLATE2(GlobalToLocalCosseratMapping, TIn, TOut),
             SOFA_TEMPLATE2(sofa::core::Mapping, TIn, TOut));

  typedef TIn In;
  typedef TOut Out;
  typedef In InDataTypes;
  typedef Out OutDataTypes;
  typedef typename InDataTypes::VecCoord InVecCoord;
  typedef typename InDataTypes::VecDeriv InVecDeriv;
  typedef typename InDataTypes::Coord InCoord;
  typedef typename InDataTypes::Deriv InDeriv;
  typedef typename InDataTypes::Real Real;
  typedef typename InDataTypes::MatrixDeriv InMatrixDeriv;

  typedef typename OutDataTypes::VecCoord OutVecCoord;
  typedef typename OutDataTypes::VecDeriv OutVecDeriv;
  typedef typename OutDataTypes::MatrixDeriv OutMatrixDeriv;
  typedef typename OutDataTypes::Coord OutCoord;
  typedef typename OutDataTypes::Deriv OutDeriv;
  typedef typename OutDataTypes::Real OutReal;

  /// @name Data Fields
  /// @{
  Data<int> d_deformationAxis;
  Data<SReal> d_max;
  Data<SReal> d_min;
  Data<SReal> d_radius;
  Data<bool> d_drawMapBeam;
  Data<sofa::type::RGBAColor> d_color;
  Data<vector<int>> d_index;
  Data<unsigned int> d_baseIndex;
  Data<vector<double>> d_curv_abs_section;
  /// @}

  struct StrainResult {
    OutCoord strain;
    Eigen::Matrix4d transformMatrix;
  };

  using Inherit1::fromModel;
  using Inherit1::toModel;

  /**
   * Calcule la déformation à partir d'un quaternion
   * @param frame Vecteur contenant la position (0-2) et le quaternion (3-6)
   * @param curvAbs Abscisse courbe
   * @param gXp Matrice précédente
   * @return Structure contenant le vecteur de déformation et la matrice de
   * transformation
   */
  StrainResult getStrainFromQuat(const InCoord &frame, const double curvAbs,
                                 const Eigen::Matrix4d &gXp);
  vector<StrainResult> m_strainResults;

public:
  void init() override;
  void apply(const sofa::core::MechanicalParams *mparams,
             Data<OutVecCoord> &out, const Data<InVecCoord> &in) override;
  void applyJ(const sofa::core::MechanicalParams *mparams,
              Data<OutVecDeriv> &out, const Data<InVecDeriv> &in) override;
  void applyJT(const sofa::core::MechanicalParams *mparams,
               Data<InVecDeriv> &outForce,
               const Data<OutVecDeriv> &inForce) override;
  void applyJT(const sofa::core::ConstraintParams *cparams,
               Data<InMatrixDeriv> &out,
               const Data<OutMatrixDeriv> &in) override;

  void computeLogarithm(const double &x, const Mat4x4 &gX, Mat4x4 &log_gX);
  void computeLogarithm(const double &x, const InCoord &gX, Mat4x4 &log_gX);

protected:
  sofa::helper::ColorMap m_colorMap;

  GlobalToLocalCosseratMapping();
  ~GlobalToLocalCosseratMapping() override = default;
};

#if !defined(SOFA_COSSERAT_CPP_GlobalToLocalCosseratMapping)
extern template class SOFA_COSSERAT_API GlobalToLocalCosseratMapping<
    sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec3Types>;
#endif

} // namespace Cosserat::mapping
