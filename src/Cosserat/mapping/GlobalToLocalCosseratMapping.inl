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
#include <Cosserat/mapping/GlobalToLocalCosseratMapping.h>

#include <sofa/core/Mapping.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/gl/template.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/helper/visual/DrawTool.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <sofa/type/Quat.h>
#include <string>

namespace Cosserat::mapping {

using sofa::core::objectmodel::BaseContext;
using sofa::defaulttype::SolidTypes;
using sofa::helper::AdvancedTimer;
using sofa::helper::WriteAccessor;
using sofa::type::RGBAColor;
template <class TIn, class TOut>
GlobalToLocalCosseratMapping<TIn, TOut>::GlobalToLocalCosseratMapping():d_deformationAxis(
    initData(&d_deformationAxis, (int)1, "deformationAxis",
      "the axis in which we want to show the deformation.\n")),
    d_max(initData(&d_max, (SReal)1.0e-2, "max",
                       "the maximum of the deformation.\n")),
    d_min(initData(&d_min, (SReal)0.0, "min",
                   "the minimum of the deformation.\n")),
    d_radius(
        initData(&d_radius, (SReal)0.05, "radius",
                 "the axis in which we want to show the deformation.\n")),
        d_color(initData(&d_color,
        sofa::type::RGBAColor(40 / 255.0, 104 / 255.0, 137 / 255.0, 0.8),
        "color", "The default beam color")),
    d_curv_abs_section(
        initData(&d_curv_abs_section, "curv_abs_section",
                 "The curvilinear abscissa of the section")) {

  }

  template <class TIn, class TOut>
  GlobalToLocalCosseratMapping<TIn, TOut>::StrainResult GlobalToLocalCosseratMapping<TIn,
  TOut>::getStrainFromQuat(const InCoord& frame,
                              const double curvAbs, const Eigen::Matrix4d& gXp)
  {
    StrainResult result;
    //result.strain.resize(3);

  // Initialisation de la matrice de transformation
  Eigen::Matrix4d gX = Eigen::Matrix4d::Zero();

  // Construction du quaternion à partir des composantes du frame
  // En python nous avons l'habitude du x,y,z,w mais ici nous avons w,x,y,z
  Eigen::Quaterniond q(frame[6], frame[3], frame[4], frame[5]); // w, x, y, z
  q.normalize();

  // Conversion du quaternion en matrice de rotation
  gX.block<3,3>(0,0) = q.toRotationMatrix();

  // Ajout de la translation
  gX.block<3,1>(0,3) = Eigen::Vector3d(frame[0], frame[1], frame[2]);
  gX(3,3) = 1.0;

  if (curvAbs <= 0.0) {
    // Cas où l'abscisse curv est nulle ou négative
    result.strain = OutCoord(0.0, 0.0, 0.0);
  } else {
    Eigen::Matrix4d gXpInv = gXp.inverse();
    Eigen::Matrix4d product = gXpInv * gX;

    // Calcul du logarithme matriciel
    //@todo : vérifier si ce calcul de log correspond à celui calculé en python
    Eigen::Matrix4d xiHat = product.log().eval() / curvAbs;

    // Extraction des composantes de déformation
    result.strain[0] = xiHat(2,1);  // κ1
    result.strain[1] = xiHat(0,2);  // κ2
    result.strain[2] = xiHat(1,0);  // κ3
  }

  result.transformMatrix = gX;
  return result;
}

  template <class TIn, class TOut>
  void GlobalToLocalCosseratMapping<TIn, TOut>::init() {

  if(fromModel.empty())
  {
    msg_error() << "input1 not found" ;
    return;
  }

  if(toModel.empty())
  {
    msg_error() << "output missing" ;
    return;
  }
  // Get initial frame state
  //auto xfromData =
  //    m_global_frames->read(sofa::core::ConstVecCoordId::position());
  //const vector<OutCoord> xfrom = xfromData->getValue();
  Inherit1::init();
}

  template <class TIn, class TOut>
  void GlobalToLocalCosseratMapping<TIn, TOut>::apply(const sofa::core::MechanicalParams *mparams,
        Data< OutVecCoord>& out, const Data< InVecCoord>& in){
      // Get initial frame state
    m_strainResults.clear();

    const auto &local_in = in.getValue();
    auto &local_out = *out.beginEdit();

    auto curvAbs = d_curv_abs_section.getValue();
    auto frames_size = local_in.size();
    local_out.resize(frames_size - 1);

    for (auto i = 0; i < frames_size; i++) {
      InCoord frame = local_in[i];
      Eigen::Matrix4d gXp = Eigen::Matrix4d::Identity();
      std::vector<double> strain;

      auto strain_struct = getStrainFromQuat(frame, curvAbs[i], gXp);
      local_out[i] = strain_struct.strain;

      m_strainResults.push_back(strain_struct);
    }
    out.endEdit();
  }

    template <class TIn, class TOut>
    void GlobalToLocalCosseratMapping<TIn, TOut>::applyJ(const sofa::core::MechanicalParams *mparams,
    Data< OutVecDeriv >& out, const Data< InVecDeriv >& in) {
      // Get initial frame state
  }
  template <class TIn, class TOut>
    void GlobalToLocalCosseratMapping<TIn, TOut>::applyJT(const sofa::core::MechanicalParams *mparams, Data< InVecDeriv >& outForce,
      const Data< OutVecDeriv >& inForce) {
      // Get initial frame state
  }

  template <class TIn, class TOut>
  void GlobalToLocalCosseratMapping<TIn, TOut>::applyJT(const sofa::core::ConstraintParams *cparams, Data< InMatrixDeriv >& out,
      const Data< OutMatrixDeriv >& in) {

  }


} // namespace Cosserat::mapping
