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

/**
 * @brief Compute logarithm using SE3's native implementation
 * 
 * @param x Scaling factor
 * @param gX Transformation matrix
 * @return Mat4x4 Logarithm of the transformation
 */
/**
 * @brief Compute logarithm using SE3's native implementation
 * 
 * @param x Scaling factor
 * @param gX Transformation matrix
 * @return Mat4x4 Logarithm of the transformation
 */
n1, TIn2, TOut>::computeLogarithm(const double &x,
                                                             const Mat4x4 &gX) -> Mat4x4
{
    if (x <= std::numeric_limits<double>::epsilon()) {
        return Mat4x4::Identity();
    }

    // Convert to Eigen format for SE3
    Eigen::Matrix4d eigen_gX;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            eigen_gX(i, j) = gX[i][j];
        }
    }
    
    // Create SE3 from matrix
    SE3d transform(eigen_gX);
    
    // Compute logarithm in the Lie algebra
    Eigen::Matrix<double, 6, 1> log_vector = transform.log();
    
    // Scale by x
    log_vector *= x;
    
    // Convert back to SE3 and then to Matrix form
    SE3d scaled_result = SE3d::exp(log_vector);
    Eigen::Matrix4d eigen_result = scaled_result.matrix();
    
    // Convert back to SOFA Mat4x4
    Mat4x4 result;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            result[i][j] = eigen_result(i, j);
        }
    }
    
    return result;
}
 Eigen::Vector3d(strain_n[0], strain_n[1], strain_n[2]);
        linear_velocity = Eigen::Vector3d(strain_n[3], strain_n[4], strain_n[5]);
    }
    
    // Scale by section length
    angular_velocity *= sub_section_length;
    linear_velocity *= sub_section_length;
    
    // Create twist vector for Lie algebra
    Eigen::Matrix<double, 6, 1> twist;
    twist << angular_velocity, linear_velocity;
    
    // Use SE3 to compute exponential map
    SE3d transform = SE3d::exp(twist);
    
    // Convert to Frame format
    g_X_n = SE3ToFrame(transform);
}
    
    return result;
}

/**
 * @brief Updates velocity state using Lie group operations
 * 
 * Implements proper velocity updates using adjoint transformations and the new 
 * Lie group functionality to propagate velocities along the beam.
 * 
 * Mathematical background:
 * For a serial kinematic chain with joint velocities k_dot, the body velocity
 * at each link can be computed recursively using:
 *   V_{i+1} = Ad_{g_{i,i+1}} * V_i + k_dot_i
 * where Ad_{g_{i,i+1}} is the adjoint of the relative transformation from link i to i+1.
 */
template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::updateVelocityState(
    const vector<Deriv1>& k_dot,
    const Vec6& base_velocity)
{
    m_nodesVelocityVectors.clear();
    m_nodeAdjointEtaVectors.clear();
    m_frameAdjointEtaVectors.clear();
    m_node_coAdjointEtaVectors.clear();
    m_frame_coAdjointEtaVectors.clear();

    // First node velocity is the base velocity
    m_nodesVelocityVectors.push_back(base_velocity);

    // Update velocities along the beam
    for (size_t i = 0; i < k_dot.size(); ++i) {
        // Store results in SOFA format
        Vec6 sofa_velocity;
        for (int j = 0; j < 6; ++j) {
            sofa_velocity[j] = new_velocity(j);
        }
        m_nodesVelocityVectors.push_back(sofa_velocity);
        
        // Store adjoint matrices for later use
        Mat6x6 sofa_adjoint;
        Mat6x6 sofa_coadjoint;
        for (int row = 0; row < 6; ++row) {
            for (int col = 0; col < 6; ++col) {
                sofa_adjoint[row][col] = adjoint(row, col);
                sofa_coadjoint[row][col] = adjoint(col, row);  // Transpose for co-adjoint
            }
        }
        
        m_nodeAdjointEtaVectors.push_back(sofa_adjoint);
        m_node_coAdjointEtaVectors.push_back(sofa_coadjoint);
    }
    
    // Calculate frame velocities
    for (size_t i = 0; i < m_indicesVectors.size(); ++i) {
        size_t beam_index = m_indicesVectors[i] - 1;
        if (beam_index < m_nodesVelocityVectors.size() - 1) {
            // Interpolate velocities for frames between beam nodes
            Vec6 node_velocity = m_nodesVelocityVectors[beam_index];
            Vec6 frame_velocity;
            
            // Transform velocity from node to frame
            transformVelocity(
                m_nodesExponentialSE3Vectors[beam_index],
                node_velocity,
                m_framesExponentialSE3Vectors[i],
                frame_velocity
            );
            
            // Store frame velocity and adjoint
            // Calculate and store adjoint matrices for frames
            SE3d frame_transform = frameToSE3(m_framesExponentialSE3Vectors[i]);
            Eigen::Matrix<double, 6, 6> frame_adjoint = frame_transform.adjoint();
            
            Mat6x6 sofa_frame_adjoint;
            Mat6x6 sofa_frame_coadjoint;
            for (int row = 0; row < 6; ++row) {
                for (int col = 0; col < 6; ++col) {
                    sofa_frame_adjoint[row][col] = frame_adjoint(row, col);
                    sofa_frame_coadjoint[row][col] = frame_adjoint(col, row);
                }
            }
            
            m_frameAdjointEtaVectors.push_back(sofa_frame_adjoint);
            m_frame_coAdjointEtaVectors.push_back(sofa_frame_coadjoint);
        }
    }
}
        Eigen::Matrix<double, 6, 1> new_velocity = adjoint * current_velocity + strain_velocity;
        
        // Store results in SOFA
            Mat6x6 sofa_frame_coadjoint;
            for (int row = 0; row < 6; ++row) {
                for (int col = 0; col < 6; ++col) {
                    sofa_frame_adjoint[row][col] = frame_adjoint(row, col);
                    sofa_frame_coadjoint[row][col] = frame_adjoint(col, row);
                }
            }
            
            m_frameAdjointEtaVectors.push_back(sofa_frame_adjoint);
            m_frame_coAdjointEtaVectors.push_back(sofa_frame_coadjoint);
        }
    }
}

/**
 * @brief Transform velocity between different coordinate frames
 * 
 * Uses SE3 adjoint to transform a velocity twist from one frame to another.
 * 
 * Mathematical background:
 * Given a twist V_a in frame A and a transformation g_ab from frame A to B,
 * the twist V_b in frame B is given by:
 *   V_b = Ad_{g_ab^{-1}} * V_a
 */
template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::transformVelocity(
    const Frame& source_frame,
    const Vec6& source_velocity,
    const Frame& target_frame,
    Vec6& target_velocity)
{
    // Convert frames to SE3
    SE3d source_transform = frameToSE3(source_frame);
    SE3d target_transform = frameToSE3(target_frame);
    
    // Compute relative transformation
    SE3d relative_transform = target_transform.inverse().compose(source_transform);
    
    // Get adjoint matrix
/**
 * @
    Eigen::Matrix<double, 6, 1> twist;
    twist << angular_velocity, linear_velocity;
    
    // Use SE3 Lie group to compute the exponential map
    SE3d transform = SE3d::exp(twist);
    
    // Convert to Frame format for SOFA
    g_X_n = SE3ToFrame(transform);
}
    // Copy to the SOFA Mat6x6 format
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            adjoint[i][j] = adjointMatrix(i, j);
        }
    }
}

/**
 * @brief Compute the co-adjoint matrix of a transformation frame
 * 
 * The co-adjoint matrix is the transpose of the adjoint matrix and is used
 * to transform wrenches (force-torque pairs) between different coordinate frames.
 * 
 * Mathematical background:
 * The co-adjoint is a 6×6 matrix that can be written in block form as:
 * 
 * Ad_g^* = [ R^T   [t]^T R^T ]
 *          [ 0     R^T      ]
 * 
 * where R is the rotation matrix, t is the translation vector, and [t] is the
 * skew-symmetric matrix formed from t.
 */
template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeCoAdjoint(const Frame &frame,
                                                             Mat6x6 &co_adjoint) {
    // Create an SE3 transformation from the Frame
    Eigen::Quaterniond quat(
        frame.getOrientation()[3],  // w
        frame.getOrientation()[0],  // x
        frame.getOrientation()[1],  // y
        frame.getOrientation()[2]   // z
    );
    
    Eigen::Vector3d trans(
        frame.getCenter()[0],
        frame.getCenter()[1], 
        frame.getCenter()[2]
    );
    
    SE3d transform(SO3d(quat), trans);
    
    // Get the co-adjoint matrix (transpose of adjoint)
    Eigen::Matrix<double, 6, 6> coAdjointMatrix = transform.adjoint().transpose();
    
    // Copy to the SOFA Mat6x6 format
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            co_adjoint[i][j] = coAdjointMatrix(i, j);
        }
    }
}

/**
 * @brief Compute the adjoint representation of a twist vector
 * 
 * This function computes the matrix representation of the adjoint action
 * of a twist vector. This matrix can be used to compute the Lie bracket
 * of two twist vectors.
 * 
 * Mathematical background:
 * For a twist ξ = (ω, v), the adjoint representation ad_ξ is a 6×6 matrix
 * that can be written in block form as:
 * 
 * ad_ξ = [ [ω]  0  ]
 *        [ [v]  [ω] ]
 * 
 * where [ω] and [v] are the skew-symmetric matrices formed from ω and v.
 */
template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeAdjoint(const Vec6 &twist,
                                                           Mat6x6 &adjoint)
{
    // Extract angular and linear velocity components
    Eigen::Vector3d omega(twist[0], twist[1], twist[2]);
    Eigen::Vector3d v(twist[3], twist[4], twist[5]);
    
    // Create the skew-symmetric matrices
    Eigen::Matrix3d omega_hat = SO3d::hat(omega);
    Eigen::Matrix3d v_hat = SO3d::hat(v);
    
    // Build the adjoint matrix
    Eigen::Matrix<double, 6, 6> ad;
    ad.setZero();
    
    // Top-left block: [ω]
    ad.block<3, 3>(0, 0) = omega_hat;
    
    // Bottom-left block: [v]
    ad.block<3, 3>(3, 0) = v_hat;
    
    // Bottom-right block: [ω]
    ad.block<3, 3>(3, 3) = omega_hat;
    
    // Copy to the SOFA Mat6x6 format
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            adjoint[i][j] = ad(i, j);
        }
    }
}
 angular and linear velocity are provided
        angular_velocity = Eigen::Vector3d(strain_n[0], strain_n[1], strain_n[2]);
        linear_velocity = Eigen::Vector3d(strain_n[3], strain_n[4], strain_n[5]);
    }
    
    // Scale by section length
    angular_velocity *= sub_section_length;
    linear_velocity *= sub_section_length;
    
    // Create twist vector (6D) for the Lie algebra element
    Eigen::Matrix<double, 6, 1> twist;
    twist << angular_velocity, linear_velocity;
    
    // Use SE3 Lie group to compute the exponential map
    SE3d transformation = SE3d::exp(twist);
    
    // Convert to Frame format for SOFA
    Eigen::Quaterniond quaternion = transformation.rotation().quaternion();
    Eigen::Vector3d translation = transformation.translation();
    
    // Update the output frame
    g_X_n = Frame(
        Vec3(translation[0], translation[1], translation[2]),
        Quat<SReal>(quaternion.w(), quaternion.x(), quaternion.y(), quaternion.z())
    );
}
    Vec3::Map(&(frame_i.getOrientation()[0])) = quaternion.coeffs();
    for (auto i = 0; i < 3; i++) {
        frame_i.getCenter()[i] = translation[i];
    }
}
            m_indicesVectorsDraw.emplace_back(sectionIndex);
        }

        // Fill the vector m_framesLengthVectors with the distance
        // between frame(output) and the closest beam node toward the base
        m_framesLengthVectors.emplace_back(
                    curv_abs_frames[i] - curv_abs_section[m_indicesVectors.back() - 1]);
    }

    for (auto j = 0; j < sz - 1; ++j)
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
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeExponentialSE3(const double &sub_section_length,
                                                                  const Coord1 &strain_n,
                                                                  Frame &g_X_n)
{
    const auto I4 = Mat4x4::Identity();

    // Get the angular part of the strain
    Vec3 k = Vec3(strain_n(0), strain_n(1), strain_n(2));
    SReal theta = k.norm();

    SE3 _g_X;
    se3 Xi_hat_n = buildXiHat(strain_n);

    //todo: change double to Real
    if (theta <= std::numeric_limits<double>::epsilon())
    {
        _g_X = I4 + sub_section_length * Xi_hat_n;
    }
    else
    {
        double scalar1 =
                (1.0 - std::cos(sub_section_length * theta)) / std::pow(theta, 2);
        double scalar2 = (sub_section_length * theta - std::sin(sub_section_length * theta)) /
                         std::pow(theta, 3);
        // Taylor expansion of exponential
        _g_X = I4 + sub_section_length * Xi_hat_n + scalar1 * Xi_hat_n * Xi_hat_n +
               scalar2 * Xi_hat_n * Xi_hat_n * Xi_hat_n;
    }

    Mat3x3 M;
    _g_X.getsub(0, 0, M); // get the rotation matrix

    // convert the rotation 3x3 matrix to a quaternion
    Quat<SReal> R; R.fromMatrix(M);
    g_X_n = Frame(Vec3(_g_X(0, 3), _g_X(1, 3), _g_X(2, 3)), R);
}

// Fill exponential vectors
template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::updateExponentialSE3(
        const vector<Coord1> &strain_state)
{
    auto curv_abs_frames = getReadAccessor(d_curv_abs_frames);

    m_framesExponentialSE3Vectors.clear();
    m_nodesExponentialSE3Vectors.clear();
    m_nodesLogarithmSE3Vectors.clear();

    const auto sz = curv_abs_frames.size();
    // Compute exponential at each frame point
    for (auto i = 0; i < sz; ++i)
    {
      Frame g_X_frame_i;

        const Coord1 strain_n =
            strain_state[m_indicesVectors[i] - 1]; // Cosserat reduce coordinates (strain)

        // the size varies from 3 to 6
        // The distance between the frame node and the closest beam node toward the base
        const SReal sub_section_length = m_framesLengthVectors[i];
        computeExponentialSE3(sub_section_length, strain_n, g_X_frame_i);
        m_framesExponentialSE3Vectors.push_back(g_X_frame_i);

        msg_info()
                << "_________________" << i << "_________________________" << msgendl
                << "x :" << sub_section_length << "; strain :" << strain_n << msgendl
                << "m_framesExponentialSE3Vectors :" << g_X_frame_i;
    }

    // Compute the exponential on the nodes
    m_nodesExponentialSE3Vectors.push_back(
        Frame(Vec3(0.0, 0.0, 0.0),
                          Quat(0., 0., 0., 1.))); // The first node.
    //todo : merge this section with the previous one
    for (unsigned int j = 0; j < strain_state.size(); ++j)
    {
        Coord1 strain_n = strain_state[j];
        const SReal section_length = m_beamLengthVectors[j];

        Frame g_X_node_j;
        computeExponentialSE3(section_length, strain_n, g_X_node_j);
        m_nodesExponentialSE3Vectors.push_back(g_X_node_j);

        msg_info()
                << "_________________Beam Node Expo___________________" << msgendl
                << "Node m_framesExponentialSE3Vectors :" << g_X_node_j << msgendl
                << "_________________Beam Node Expo___________________";

    }
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeAdjoint(const Frame &frame,
                                                           TangentTransform &adjoint)
{
    Mat3x3 R = extractRotMatrix(frame);
    Vec3 u = frame.getOrigin();
    Mat3x3 tilde_u = getTildeMatrix(u);
    Mat3x3 tilde_u_R = tilde_u * R;
    buildAdjoint(R, tilde_u_R, adjoint);
}

template <class TIn1, class TIn2, class TOut>
void BaseCosseratMapping<TIn1, TIn2, TOut>::computeCoAdjoint(const Frame &frame,
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
    auto d_x1 = m_strain_state->read(sofa::core::vec_id::read_access::position);

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

    Frame outTransform;
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
Mat3x3 BaseCosseratMapping<TIn1, TIn2, TOut>::extractRotMatrix(const Frame &frame) {

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
auto BaseCosseratMapping<TIn1, TIn2, TOut>::buildProjector(const Frame &T)
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
        const Frame &T) -> Mat4x4 {
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
