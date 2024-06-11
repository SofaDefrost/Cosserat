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
#define SOFA_COSSERAT_CPP_BaseCosseratMapping
#include <Cosserat/mapping/BaseCosseratMapping.inl>

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

namespace Cosserat::mapping
{
using namespace sofa::defaulttype;

template <>
BaseCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>::se3
BaseCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>::buildXiHat(const Coord1 &strain_i)
{
  se3 Xi;

  Xi[0][1] = -strain_i(2);
  Xi[0][2] = strain_i[1];
  Xi[1][2] = -strain_i[0];

  Xi[1][0] = -Xi(0, 1);
  Xi[2][0] = -Xi(0, 2);
  Xi[2][1] = -Xi(1, 2);

  Xi[0][3] = 1.0;

  for (unsigned int i = 0; i < 3; i++)
    Xi[i][3] += strain_i(i + 3);

  //    se3 = [
  //        0               -screw(3)   screw(2)        screw(4);
  //        screw(3)        0           -screw(1)       screw(5);
  //        -screw(2)   screw(1)        0               screw(6);
  //        0                   0                 0                 0];

  return Xi;
}

template <>
void BaseCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>::computeExponentialSE3(
    const double &curv_abs_x_n, const Coord1 &strain_n, Transform &Trans) {
  Matrix4 I4;
  I4.identity();

  // Get the angular part of the
  Vec3 k = Vec3(strain_n(0), strain_n(1), strain_n(2));
  SReal theta = k.norm(); //

  SE3 g_X_n;
  se3 Xi_hat_n = buildXiHat(strain_n);

  msg_info() << "matrix Xi : " << Xi_hat_n;

  if (theta <= std::numeric_limits<double>::epsilon()) {
    g_X_n = I4 + curv_abs_x_n * Xi_hat_n;
  } else {
    double scalar1 =
        (1.0 - std::cos(curv_abs_x_n * theta)) / std::pow(theta, 2);
    double scalar2 = (curv_abs_x_n * theta - std::sin(curv_abs_x_n * theta)) /
                     std::pow(theta, 3);
    g_X_n = I4 + curv_abs_x_n * Xi_hat_n + scalar1 * Xi_hat_n * Xi_hat_n +
            scalar2 * Xi_hat_n * Xi_hat_n * Xi_hat_n;
  }

  msg_info() << "matrix g_X : " << g_X_n;

  sofa::type::Mat3x3 M;
  g_X_n.getsub(0, 0, M); // get the rotation matrix

  // convert the rotation 3x3 matrix to a quaternion
  sofa::type::Quat<Real> R;
  R.fromMatrix(M);
  Trans = Transform(Vec3(g_X_n(0, 3), g_X_n(1, 3), g_X_n(2, 3)), R);
}

template <>
void BaseCosseratMapping<Vec6Types, Rigid3Types, Rigid3Types>::computeTangExp(
    double &curv_abs_n, const Coord1 &strain_i, Mat6x6 &TgX) {

  SReal theta = Vec3(strain_i(0), strain_i(1), strain_i(2))
                    .norm(); // Sometimes this is computed over all strain
  Matrix3 tilde_k =
      getTildeMatrix(Vec3(strain_i(0), strain_i(1), strain_i(2)));

  /* Younes @23-11-27
  old version
  @Todo ???? is p the linear deformation? If so, why didn't I just put a zero
  vector in place of p and the first element of p is equal to 1? Matrix3 tilde_p
  = getTildeMatrix(type::Vec3(1.0, 0.0, 0.0)); Using the new version does not
  bring any difference in my three reference scenes, but need more investogation
  #TECHNICAL_DEBT
  */
  // TODO(dmarchal: 2024/06/07) could the debt by solved ?
  Matrix3 tilde_q =
      getTildeMatrix(Vec3(strain_i(3), strain_i(4), strain_i(5)));

  Mat6x6 ad_Xi;
  buildAdjoint(tilde_k, tilde_q, ad_Xi);

  Mat6x6 Id6 = Mat6x6::Identity();

  if (theta <= std::numeric_limits<double>::epsilon()) {
    double scalar0 = std::pow(curv_abs_n, 2) / 2.0;
    TgX = curv_abs_n * Id6 + scalar0 * ad_Xi;
  } else {
    double scalar1 = (4.0 - 4.0 * cos(curv_abs_n * theta) -
                      curv_abs_n * theta * sin(curv_abs_n * theta)) /
                     (2.0 * theta * theta);
    double scalar2 = (4.0 * curv_abs_n * theta +
                      curv_abs_n * theta * cos(curv_abs_n * theta) -
                      5.0 * sin(curv_abs_n * theta)) /
                     (2.0 * theta * theta * theta);
    double scalar3 = (2.0 - 2.0 * cos(curv_abs_n * theta) -
                      curv_abs_n * theta * sin(curv_abs_n * theta)) /
                     (2.0 * theta * theta * theta * theta);
    double scalar4 = (2.0 * curv_abs_n * theta +
                      curv_abs_n * theta * cos(curv_abs_n * theta) -
                      3.0 * sin(curv_abs_n * theta)) /
                     (2.0 * theta * theta * theta * theta * theta);

    TgX = curv_abs_n * Id6 + scalar1 * ad_Xi + scalar2 * ad_Xi * ad_Xi +
          scalar3 * ad_Xi * ad_Xi * ad_Xi +
          scalar4 * ad_Xi * ad_Xi * ad_Xi * ad_Xi;
  }
}


// Register in the Factory
int BaseCosseratMappingClass =
    sofa::core::RegisterObject(
        "Set the positions and velocities of points attached to a rigid parent")
        .add<BaseCosseratMapping<sofa::defaulttype::Vec3Types,
                          sofa::defaulttype::Rigid3Types,
                          sofa::defaulttype::Rigid3Types>>()
        .add<BaseCosseratMapping<sofa::defaulttype::Vec6Types,
                          sofa::defaulttype::Rigid3Types,
                          sofa::defaulttype::Rigid3Types>>();

template class SOFA_COSSERAT_API
    BaseCosseratMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
                 sofa::defaulttype::Rigid3Types>;
template class SOFA_COSSERAT_API
    BaseCosseratMapping<sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types,
                 sofa::defaulttype::Rigid3Types>;

} // namespace cosserat::mapping
