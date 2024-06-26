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

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/Quat.h>

#include <string>

// To go further =>
// https://www.mathworks.com/matlabcentral/fileexchange/83038-sorosim/

namespace Cosserat::mapping
{

using sofa::helper::getReadAccessor;
using sofa::type::Vec6;
using sofa::type::Vec3;
using sofa::type::Quat;

template <class TIn1, class TIn2, class TOut>
BaseCosseratMapping<TIn1, TIn2, TOut>::BaseCosseratMapping()
    // TODO(dmarchal: 2024/06/12): please add the help comments !
    : d_curv_abs_section(initData(&d_curv_abs_section, "curv_abs_input",
                                  " need to be com....")),
      d_curv_abs_frames(initData(&d_curv_abs_frames, "curv_abs_output",
                                 " need to be com....")),
      d_debug(initData(&d_debug, false, "debug", "printf for the debug")),
      m_indexInput(0) {}


template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::init()
{
    m_fromModel1 = nullptr;
    m_fromModel2 = nullptr;
    m_toModel = nullptr;

    if(fromModels1.empty())
    {
        msg_error() << "input1 not found" ;
        return;
    }

    if(fromModels2.empty())
    {
        msg_error() << "input2 not found" ;
        return;
    }

    if(toModels.empty())
    {
        msg_error() << "output missing" ;
        return;
    }

    m_fromModel1 = fromModels1[0];
    m_fromModel2 = fromModels2[0];
    m_toModel = toModels[0];

    // Fill the initial vector
    auto xfromData = m_toModel->read(sofa::core::ConstVecCoordId::position());
    const vector<OutCoord> xfrom = xfromData->getValue();

    m_vecTransform.clear();
    for (unsigned int i = 0; i < xfrom.size(); i++)
        m_vecTransform.push_back(xfrom[i]);

    initializeFrames();
    doBaseCosseratInit();
    Inherit1::init();
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::initializeFrames()
{
    // For each frame in the global frame, find the segment of the beam to which
    // it is attached. Here we only use the information from the curvilinear
    // abscissa of each frame.
    auto curv_abs_section = getReadAccessor(d_curv_abs_section);
    auto curv_abs_frames = getReadAccessor(d_curv_abs_frames);

    msg_info()
            << " curv_abs_section " << curv_abs_section.size() << "; curv_abs_frames: " << curv_abs_frames.size();

    m_indicesVectors.clear();
    m_framesLengthVectors.clear();
    m_beamLengthVectors.clear();
    m_indicesVectorsDraw.clear();

    size_t sz = curv_abs_frames.size();
    size_t input_index = 1;
    for (size_t i = 0; i < sz; ++i)
    {
        if (curv_abs_section[input_index] > curv_abs_frames[i])
        {
            m_indicesVectors.emplace_back(input_index);
            m_indicesVectorsDraw.emplace_back(input_index); // maybe I shouldn't do this here !!!
        } else if (curv_abs_section[input_index] == curv_abs_frames[i])
        {
            m_indicesVectors.emplace_back(input_index);
            input_index++;
            m_indicesVectorsDraw.emplace_back(input_index);
        } else {
            input_index++;
            m_indicesVectors.emplace_back(input_index);
            m_indicesVectorsDraw.emplace_back(input_index);
        }

        // Fill the vector m_framesLengthVectors with the distance
        // between frame(output) and the closest beam node toward the base
        m_framesLengthVectors.emplace_back(
                    curv_abs_frames[i] - curv_abs_section[m_indicesVectors.back() - 1]);
    }

    for (size_t j = 0; j < sz - 1; ++j)
    {
        m_beamLengthVectors.emplace_back(curv_abs_section[j + 1] - curv_abs_section[j]);
    }

    msg_info()
            << "m_indicesVectors : " << m_indicesVectors << msgendl
            << "m_framesLengthVectors : " << msgendl
            << "m_BeamLengthVectors : " << msgendl;
}

auto buildXiHat(const Vec3& strain_i) -> se3
{
    se3 Xi_hat;

    Xi_hat[0][1] = -strain_i[2];
    Xi_hat[0][2] =  strain_i[1];
    Xi_hat[1][2] = -strain_i[0];

    Xi_hat[1][0] = -Xi_hat(0, 1);
    Xi_hat[2][0] = -Xi_hat(0, 2);
    Xi_hat[2][1] = -Xi_hat(1, 2);

    //@TODO:  Why this , if q = 0 ????
    Xi_hat[0][3] = 1.0;
    return Xi_hat;
}

auto buildXiHat(const Vec6& strain_i) -> se3
{
    se3 Xi = buildXiHat(Vec3(strain_i(0), strain_i(1), strain_i(2)));

    for (unsigned int i = 0; i < 3; i++)
      Xi[i][3] += strain_i(i + 3);

    return Xi;
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeExponentialSE3(const double &curv_abs_x_n,
                                                                  const Coord1 &strain_n,
                                                                  Transform &g_X_n)
{
    const auto I4 = Mat4x4::Identity();

    // Get the angular part of the strain
    Vec3 k = Vec3(strain_n(0), strain_n(1), strain_n(2));
    SReal theta = k.norm();

    SE3 _g_X;
    se3 Xi_hat_n = buildXiHat(strain_n);

    msg_info() << "matrix Xi : " << Xi_hat_n;

    if (theta <= std::numeric_limits<double>::epsilon()) {
        _g_X = I4 + curv_abs_x_n * Xi_hat_n;
    } else {
        double scalar1 =
                (1.0 - std::cos(curv_abs_x_n * theta)) / std::pow(theta, 2);
        double scalar2 = (curv_abs_x_n * theta - std::sin(curv_abs_x_n * theta)) /
                         std::pow(theta, 3);
        _g_X = I4 + curv_abs_x_n * Xi_hat_n + scalar1 * Xi_hat_n * Xi_hat_n +
               scalar2 * Xi_hat_n * Xi_hat_n * Xi_hat_n;
    }

    msg_info() << "matrix _g_X : " << _g_X;

    Mat3x3 M;
    _g_X.getsub(0, 0, M); // get the rotation matrix

    // convert the rotation 3x3 matrix to a quaternion
    Quat<SReal> R;
    R.fromMatrix(M);
    g_X_n = Transform(Vec3(_g_X(0, 3), _g_X(1, 3), _g_X(2, 3)), R);
}

// Fill exponential vectors
template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::updateExponentialSE3(
        const vector<Coord1> &inDeform)
{
    auto curv_abs_frames = getReadAccessor(d_curv_abs_frames);

    m_framesExponentialSE3Vectors.clear();
    m_nodesExponentialSE3Vectors.clear();
    m_nodesLogarithmeSE3Vectors.clear();

    const unsigned int sz = curv_abs_frames.size();

    // Compute exponential at each frame point
    for (size_t i = 0; i < sz; ++i)
    {
        Transform g_X_frame_i;

        const Coord1 strain_n = inDeform[m_indicesVectors[i] - 1]; // Cosserat reduce coordinates (strain)

        // the size varies from 1 to 6
        // The distance between the frame and the closest beam node toward the base
        const SReal curv_abs_x = m_framesLengthVectors[i];
        computeExponentialSE3(curv_abs_x, strain_n, g_X_frame_i);
        m_framesExponentialSE3Vectors.push_back(g_X_frame_i);

        msg_info()
                << "_________________" << i << "_________________________" << msgendl
                << "x :" << curv_abs_x << "; strain :" << strain_n << msgendl
                << "m_framesExponentialSE3Vectors :" << g_X_frame_i;
    }

    // Compute the exponential on the nodes
    m_nodesExponentialSE3Vectors.push_back(
                Transform(Vec3(0.0, 0.0, 0.0),
                          Quat(0., 0., 0., 1.))); // The first node.

    for (unsigned int j = 0; j < inDeform.size(); ++j)
    {
        Coord1 strain_n = inDeform[j];
        const SReal curv_abs_x = m_beamLengthVectors[j];

        Transform g_X_node_j;
        computeExponentialSE3(curv_abs_x, strain_n, g_X_node_j);
        m_nodesExponentialSE3Vectors.push_back(g_X_node_j);

        msg_info()
                << "_________________Beam Node Expo___________________" << msgendl
                << "Node m_framesExponentialSE3Vectors :" << g_X_node_j << msgendl
                << "_________________Beam Node Expo___________________";

    }
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeAdjoint(const Transform &frame,
                                                           TangentTransform &adjoint)
{
    Mat3x3 R = extractRotMatrix(frame);
    Vec3 u = frame.getOrigin();
    Mat3x3 tilde_u = getTildeMatrix(u);
    Mat3x3 tilde_u_R = tilde_u * R;
    buildAdjoint(R, tilde_u_R, adjoint);
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeCoAdjoint(const Transform &frame,
                                                             Mat6x6 &co_adjoint) {
    Mat3x3 R = extractRotMatrix(frame);
    Vec3 u = frame.getOrigin();
    Mat3x3 tilde_u = getTildeMatrix(u);
    Mat3x3 tilde_u_R = tilde_u * R;
    buildCoAdjoint(R, tilde_u_R, co_adjoint);
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeAdjoint(const Vec6 &eta,
                                                           Mat6x6 &adjoint)
{
    Mat3x3 tildeMat1 = getTildeMatrix(Vec3(eta[0], eta[1], eta[2]));
    Mat3x3 tildeMat2 = getTildeMatrix(Vec3(eta[3], eta[4], eta[5]));
    buildAdjoint(tildeMat1, tildeMat2, adjoint);
}


template <class TIn1, class TIn2, class TOut>
auto BaseCosseratMapping<TIn1, TIn2, TOut>::computeLogarithm(const double &x,
                                                             const Mat4x4 &gX) -> Mat4x4
{
    // Compute theta before everything
    const double theta = computeTheta(x, gX);
    Mat4x4 I4 = Mat4x4::Identity();
    Mat4x4 log_gX;

    double csc_theta = 1.0 / (sin(x * theta / 2.0));
    double sec_theta = 1.0 / (cos(x * theta / 2.0));
    double cst = (1.0 / 8) * (csc_theta * csc_theta * csc_theta) * sec_theta;
    double x_theta = x * theta;
    double cos_2XTheta = cos(2.0 * x_theta);
    double cos_XTheta = cos(x_theta);
    double sin_2XTheta = sin(2.0 * x_theta);
    double sin_XTheta = sin(x_theta);

    if (theta <= std::numeric_limits<double>::epsilon())
        log_gX = I4;
    else {
        log_gX = cst * ((x_theta * cos_2XTheta - sin_XTheta) * I4 -
                        (x_theta * cos_XTheta + 2.0 * x_theta * cos_2XTheta -
                         sin_XTheta - sin_2XTheta) *
                        gX +
                        (2.0 * x_theta * cos_XTheta + x_theta * cos_2XTheta -
                         sin_XTheta - sin_2XTheta) *
                        (gX * gX) -
                        (x_theta * cos_XTheta - sin_XTheta) * (gX * gX * gX));
    }

    return log_gX;
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::updateTangExpSE3(
        const vector<Coord1> &inDeform) {

    // Curv abscissa of nodes and frames
    auto curv_abs_section = getReadAccessor(d_curv_abs_section);
    auto curv_abs_frames = getReadAccessor(d_curv_abs_frames);

    unsigned int sz = curv_abs_frames.size();
    m_framesTangExpVectors.resize(sz);

    // Compute tangExpo at frame points
    for (unsigned int i = 0; i < sz; i++)
    {
        TangentTransform temp;

        Coord1 strain_frame_i = inDeform[m_indicesVectors[i] - 1];
        double curv_abs_x_i = m_framesLengthVectors[i];
        computeTangExp(curv_abs_x_i, strain_frame_i, temp);

        m_framesTangExpVectors[i] = temp;

        msg_info()
                <<  "x :" << curv_abs_x_i << "; k :" << strain_frame_i << msgendl
                 << "m_framesTangExpVectors :" << m_framesTangExpVectors[i];
    }

    // Compute the TangExpSE3 at the nodes
    m_nodesTangExpVectors.clear();
    TangentTransform tangExpO;
    tangExpO.clear();
    m_nodesTangExpVectors.push_back(tangExpO);

    for (size_t j = 1; j < curv_abs_section.size(); j++) {
        Coord1 strain_node_i = inDeform[j - 1];
        double x = m_beamLengthVectors[j - 1];
        TangentTransform temp;
        temp.clear();
        computeTangExp(x, strain_node_i, temp);
        m_nodesTangExpVectors.push_back(temp);
    }
    msg_info() << "Node TangExpo : " << m_nodesTangExpVectors;
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeTangExp(double &curv_abs_n,
                                                           const Coord1 &strain_i,
                                                           Mat6x6 &TgX)
{
    if constexpr( Coord1::static_size == 3 )
        computeTangExpImplementation(curv_abs_n, Vec6(strain_i(0),strain_i(1),strain_i(2),0,0,0), TgX);
    else
        computeTangExpImplementation(curv_abs_n, strain_i, TgX);
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeTangExpImplementation(double &curv_abs_n,
                                                                         const Vec6 &strain_i,
                                                                         Mat6x6 &TgX)
{
    SReal theta = Vec3(strain_i(0), strain_i(1), strain_i(2)).norm();
    Mat3x3 tilde_k = getTildeMatrix(Vec3(strain_i(0), strain_i(1), strain_i(2)));
    Mat3x3 tilde_q = getTildeMatrix(Vec3(strain_i(3), strain_i(4), strain_i(5)));

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

template <class TIn1, class TIn2, class TOut>
[[maybe_unused]] Vec6
BaseCosseratMapping<TIn1, TIn2, TOut>::computeETA(const Vec6 &baseEta,
                                                  const vector<Deriv1> &k_dot,
                                                  const double abs_input)
{

    // Get the positions from model 0. This function returns the position wrapped in a Data<>
    auto d_x1 = m_fromModel1->read(sofa::core::ConstVecCoordId::position());

    // To access the actual content (in this case position) from a data, we have to use
    // a read accessor that insures the data is updated according to DDGNode state
    auto x1 = getReadAccessor(*d_x1);

    // Same as for x1, query a read accessor so we can access the content of d_curv_abs_section
    auto curv_abs_input = getReadAccessor(d_curv_abs_section);

    auto& kdot = k_dot[m_indexInput];
    Vec6 Xi_dot {kdot[0], kdot[1], kdot[2],
                 0,0,0};

    // if m_indexInput is == 0
    double diff0 = abs_input;
    double _diff0 = -abs_input;

    if (m_indexInput != 0)
    {
        diff0 = abs_input - curv_abs_input[m_indexInput - 1];
        _diff0 = curv_abs_input[m_indexInput - 1] - abs_input;
    }

    Transform outTransform;
    computeExponentialSE3(_diff0, x1[m_indexInput], outTransform);

    TangentTransform adjointMatrix;
    computeAdjoint(outTransform, adjointMatrix);

    TangentTransform tangentMatrix;
    computeTangExp(diff0, x1[m_indexInput], tangentMatrix);

    return adjointMatrix * (baseEta + tangentMatrix * Xi_dot);
}


template <class TIn1, class TIn2, class TOut>
double BaseCosseratMapping<TIn1, TIn2, TOut>::computeTheta(const double &x,
                                                           const Mat4x4 &gX) {
    double Tr_gx = sofa::type::trace(gX);

    if (x > std::numeric_limits<double>::epsilon())
        return (1.0 / x) * std::acos((Tr_gx / 2.0) - 1);

    return 0.0;
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::printMatrix(const Mat6x6 R) {
    // TODO(dmarchal: 2024/06/07): Remove the use of printf in addition to
    // reconsider the implementation of common utility functions in instance
    // method.
    for (unsigned int k = 0; k < 6; k++) {
        for (unsigned int i = 0; i < 6; i++)
            printf("  %lf", R[k][i]);
        printf("\n");
    }
}

template <class TIn1, class TIn2, class TOut>
Mat3x3 BaseCosseratMapping<TIn1, TIn2, TOut>::extractRotMatrix(const Transform &frame) {

    Quat q = frame.getOrientation();

    // TODO(dmarchal: 2024/06/07) The following code should probably become
    // utility function building a 3x3 matix from a quaternion should probably
    // does not need this amount of code.
    SReal R[4][4];
    q.buildRotationMatrix(R);
    Mat3x3 mat;
    for (unsigned int k = 0; k < 3; k++)
        for (unsigned int i = 0; i < 3; i++)
            mat[k][i] = R[k][i];
    return mat;
}

template <class TIn1, class TIn2, class TOut>
auto BaseCosseratMapping<TIn1, TIn2, TOut>::buildProjector(const Transform &T)
-> TangentTransform {
    TangentTransform P;

    // TODO(dmarchal: 2024/06/07) The following code should probably become
    // utility function building a 3x3 matix from a quaternion should probably
    // does not need this amount of code.
    SReal R[4][4];
    (T.getOrientation()).buildRotationMatrix(R);
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 3; j++) {
            P[i][j + 3] = R[i][j];
            P[i + 3][j] = R[i][j];
        }
    }
    return P;
}

template <class TIn1, class TIn2, class TOut>
auto BaseCosseratMapping<TIn1, TIn2, TOut>::getTildeMatrix(const Vec3 &u)
-> Mat3x3 {
    sofa::type::Matrix3 tild;
    tild[0][1] = -u[2];
    tild[0][2] = u[1];
    tild[1][2] = -u[0];

    tild[1][0] = -tild[0][1];
    tild[2][0] = -tild[0][2];
    tild[2][1] = -tild[1][2];
    return tild;
}

template <class TIn1, class TIn2, class TOut>
auto BaseCosseratMapping<TIn1, TIn2, TOut>::buildAdjoint(const Mat3x3 &A,
                                                         const Mat3x3 &B,
                                                         Mat6x6 &Adjoint) -> void {
    Adjoint.clear();
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Adjoint[i][j] = A[i][j];
            Adjoint[i + 3][j + 3] = A[i][j];
            Adjoint[i + 3][j] = B[i][j];
        }
    }
}

template <class TIn1, class TIn2, class TOut>
auto BaseCosseratMapping<TIn1, TIn2, TOut>::buildCoAdjoint(const Mat3x3 &A,
                                                           const Mat3x3 &B,
                                                           Mat6x6 &coAdjoint) -> void {
    coAdjoint.clear();
    for (unsigned int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            coAdjoint[i][j] = A[i][j];
            coAdjoint[i + 3][j + 3] = A[i][j];

            // TODO(dmarchal: 2024/06/07) if co-adjoint is juste this single change
            // (the j+3)
            coAdjoint[i][j + 3] = B[i][j];
        }
    }
}

template <class TIn1, class TIn2, class TOut>
auto BaseCosseratMapping<TIn1, TIn2, TOut>::convertTransformToMatrix4x4(
        const Transform &T) -> Mat4x4 {
    Mat4x4 M = Mat4x4::Identity();
    Mat3x3 R = extractRotMatrix(T);
    Vec3 trans = T.getOrigin();

    for (unsigned int i = 0; i < 3; i++)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            M(i, j) = R[i][j];
            M(i, 3) = trans[i];
        }
    }
    return M;
}

template <class TIn1, class TIn2, class TOut>
auto BaseCosseratMapping<TIn1, TIn2, TOut>::piecewiseLogmap(const _SE3 &g_x) -> Vec6 {
    _SE3 Xi_hat;

    double x = 1.0;
    double theta = std::acos(g_x.trace() / 2.0 - 1.0);

    if (theta == 0) {
        Xi_hat = 1.0 / x * (g_x - Matrix4d::Identity());
    } else {
        double x_theta = x * theta;
        double sin_x_theta = std::sin(x_theta);
        double cos_x_theta = std::cos(x_theta);
        double t3 = 2 * sin_x_theta * cos_x_theta;
        double t4 = 1 - 2 * sin_x_theta * sin_x_theta;
        double t5 = x_theta * t4;

        Matrix4d gp2 = g_x * g_x;
        Matrix4d gp3 = gp2 * g_x;

        Xi_hat = 1.0 / x *
                 (0.125 *
                  (1.0 / std::sin(x_theta / 2.0) / std::sin(x_theta / 2.0) /
                   std::sin(x_theta / 2.0)) *
                  std::cos(x_theta / 2.0) *
                  ((t5 - sin_x_theta) * Matrix4d::Identity() -
                   (x_theta * cos_x_theta + 2 * t5 - sin_x_theta - t3) * g_x +
                   (2 * x_theta * cos_x_theta + t5 - sin_x_theta - t3) * gp2 -
                   (x_theta * cos_x_theta - sin_x_theta) * gp3));
    }

    Vec6 xci = Vec6(Xi_hat(2, 1), Xi_hat(0, 2), Xi_hat(1, 0), Xi_hat(0, 3),
                    Xi_hat(1, 3), Xi_hat(2, 3));
    return xci;
}

} // namespace cosserat::mapping
