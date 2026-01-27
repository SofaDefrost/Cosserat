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
#include <Cosserat/types.h>
#include <Cosserat/liegroups/SO3.h>
#include <Cosserat/liegroups/SE3.h>
#include <Cosserat/liegroups/Utils.h>

#include <sofa/core/Multi2Mapping.h>

namespace Cosserat::mapping
{

// Use a private namespace so we are not polluting the Cosserat::mapping.
namespace
{
using namespace std;
using namespace Eigen;
using sofa::defaulttype::SolidTypes;
using sofa::type::Mat6x6;
using sofa::type::Mat4x4;
using sofa::type::Mat3x3;

using std::get;
using sofa::type::vector;
using sofa::type::Vec3;
using sofa::type::Vec6;
using sofa::type::Mat;

// Modern Lie group implementations
using SO3d = Cosserat::SO3<double>;
using SE3d = Cosserat::SE3<double>;

// Legacy types for backward compatibility
using SE3_Legacy = sofa::type::Matrix4; ///< Legacy Matrix4 representation
using se3_Legacy = sofa::type::Matrix4; ///< Legacy Matrix4 representation
using _se3 = Eigen::Matrix4d;
using _SE3 = Eigen::Matrix4d;

using Cosserat::type::Frame;
using Cosserat::type::TangentTransform;
using Cosserat::type::RotMat;

}
}
/*!
 * \class BaseCosseratMapping
 * @brief Base class for Cosserat rod mappings in SOFA framework
 *
 * This class provides the foundation for implementing Cosserat rod mappings,
 * which are used to map between different representations of a Cosserat rod's
 * configuration and deformation.
 *
 * @tparam TIn1 The first input type for the mapping
 * @tparam TIn2 The second input type for the mapping
 * @tparam TOut The output type for the mapping
 */
template <class TIn1, class TIn2, class TOut>
class BaseCosseratMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut>
{
public:
    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE3(BaseCosseratMapping, TIn1, TIn2, TOut),
                        SOFA_TEMPLATE3(sofa::core::Multi2Mapping,TIn1, TIn2, TOut));

    typedef TIn1 In1;
    typedef TIn2 In2;
    typedef TOut Out;

    using Coord1 = sofa::Coord_t<In1>;
    using Deriv1 = sofa::Deriv_t<In1>;
    using OutCoord = sofa::Coord_t<Out>;

    /*===========COSSERAT VECTORS ======================*/
    unsigned int m_indexInput;
    vector<OutCoord> m_vecTransform;

    vector<Frame> m_framesExponentialSE3Vectors;
    vector<Frame> m_nodesExponentialSE3Vectors;
    vector<Mat4x4> m_nodesLogarithmSE3Vectors;

    vector<unsigned int> m_indicesVectors;
    vector<unsigned int> m_indicesVectorsDraw;

    vector<double> m_beamLengthVectors;
    vector<double> m_framesLengthVectors;

    vector<Vec6> m_nodesVelocityVectors;
    vector<Mat6x6> m_nodesTangExpVectors;
    vector<Mat6x6> m_framesTangExpVectors;
    vector<Vec6> m_totalBeamForceVectors;

    vector<Mat6x6> m_nodeAdjointVectors;

    /**
     * @brief Compute the adjoint representation for a transformation frame
     * 
     * The adjoint representation is used to transform twists between different
     * coordinate frames. It is also used for updating velocities.
     * 
     * @param frame The transformation frame
     * @param adjoint Output adjoint matrix as TangentTransform
     */
    /**
     * @brief Compute the co-adjoint matrix of a transformation frame
     * 
     * The co-adjoint matrix is the transpose of the adjoint matrix, used
     * to transform wrenches (force-torque pairs) between coordinate frames.
     * 
     * @param frame The transformation frame
     * @param coAdjoint Output co-adjoint matrix
     */
    void computeCoAdjoint(const Frame &frame, Mat6x6 &coAdjoint);

    /**
     * @brief Update exponential vectors for all frames and nodes
     * 
     * @param strain_state Vector of strain states
     */
    void updateExponentialSE3(const vector<Coord1> &strain_state);
    
    /**
     * @brief Update tangent exponential vectors
     * 
     * @param inDeform Vector of deformations
     */
    void updateTangExpSE3(const vector<Coord1> &inDeform);

    /**
     * @brief Compute tangent exponential map
     * 
     * @param x Parameter for tangent map
     * @param k Strain vector
     * @param TgX Output tangent matrix
     */
    void computeTangExp(double &x, const Coord1 &k, Mat6x6 &TgX);
    
    /**
     * @brief Implementation of tangent exponential map
     * 
     * @param x Parameter for tangent map
     * @param k Strain vector
     * @param TgX Output tangent matrix
     */
    void computeTangExpImplementation(double &x, const Vec6 &k, Mat6x6 &TgX);

    /**
     * @brief Compute eta vector for a given input
     * 
     * @param baseEta Base eta vector
     * @param k_dot Vector of strain rates
     * @param abs_input Position along the rod
     * @return Vec6 Computed eta vector
     */
    [[maybe_unused]] Vec6
    computeETA(const Vec6 &baseEta, const vector<Deriv1> &k_dot, double abs_input);
    
    /**
     * @brief Compute logarithm map using SE3
     * 
     * @param x Scaling factor
     * @param gX Transformation matrix
     * @return Mat4x4 Logarithm of the transformation
     */
    Mat4x4 computeLogarithm(const double &x, const Mat4x4 &gX);
    void computeAdjoint(const Vec6 &twist, Mat6x6 &adjoint);
    
    /**
     * @brief Updates velocity state using Lie group operations
     * 
     * Implements proper velocity updates using adjoint transformations and
     * the new Lie group functionality to propagate velocities along the beam.
     * 
     * @param k_dot Strain rates (angular and linear velocity derivatives)
     * @param base_velocity Base node velocity in body coordinates
     */
    void updateVelocityState(const vector<Deriv1>& k_dot, const Vec6& base_velocity);
    
    /**
     * @brief Transform velocity between different coordinate frames
     * 
     * Uses SE3 adjoint to transform a velocity twist from one frame to another.
     * 
     * @param source_frame Source coordinate frame
     * @param source_velocity Velocity in source frame
     * @param target_frame Target coordinate frame
     * @param target_velocity Output: velocity expressed in target frame
     */
    void transformVelocity(
        const Frame& source_frame,
        const Vec6& source_velocity,
        const Frame& target_frame,
        Vec6& target_velocity);

    /**
     * @brief Compute the angle parameter for logarithm calculation
     *
     * @param x Scaling factor
     * @param gX Transformation matrix
     * @return double The angle parameter
     */
    double computeTheta(const double &x, const Mat4x4 &gX);

    /**
     * @brief Extract rotation matrix from a Frame using SO3
     * 
     * @param frame The input Frame containing orientation as a quaternion
     * @return Mat3x3 The 3x3 rotation matrix
     */
    Mat3x3 extractRotMatrix(const Frame &frame);

    /**
     * @brief Build a projector matrix from a Frame
     *
     * @param T The transformation frame
     * @return TangentTransform The projector matrix
     */
    TangentTransform buildProjector(const Frame &T);

    /**
     * @brief Create a skew-symmetric matrix from a vector using SO3::hat
     * 
     * @param u The input 3D vector
     * @return sofa::type::Matrix3 The skew-symmetric matrix
     */
    sofa::type::Matrix3 getTildeMatrix(const Vec3 &u);

    /**
     * @brief Print a matrix using modern logging
     * 
     * Uses SOFA's message system instead of printf for better integration
     * with the logging framework.
     * 
     * @param matrix The 6x6 matrix to print
     */
    void printMatrix(const Mat6x6& matrix);
    
    /**
     * @brief Convert a SOFA Frame to an SE3 transformation
     * 
     * @param frame The input SOFA Frame
     * @return SE3d The SE3 transformation
     */
    SE3d frameToSE3(const Frame &frame);
    
    /**
     * @brief Convert an SE3 transformation to a SOFA Frame
     * 
     * @param transform The input SE3 transformation
     * @return Frame The SOFA Frame
     */
    Frame SE3ToFrame(const SE3d &transform);

    Mat4x4 convertTransformToMatrix4x4(const Frame &T);
    Vec6 piecewiseLogmap(const _SE3 &g_x);

    /*!
     * @brief Computes the rotation matrix around the X-axis using SO3
     *
     * Uses the Lie group SO3 implementation for better numerical stability
     * and consistency with other transformations.
     *
     * @param angle The rotation angle in radians
     * @return RotMat A 3x3 rotation matrix representing the rotation around the X-axis
     */
    RotMat rotationMatrixX(double angle) {
        // Create rotation using the exponential map with axis (1,0,0)
        Eigen::Vector3d axis = Eigen::Vector3d::UnitX();
        return SO3d::exp(angle * axis).matrix();
    }

    /*!
     * @brief Computes the rotation matrix around the Y-axis using SO3
     *
     * Uses the Lie group SO3 implementation for better numerical stability
     * and consistency with other transformations.
     *
     * @param angle The rotation angle in radians
     * @return RotMat A 3x3 rotation matrix representing the rotation around the Y-axis
     */
    RotMat rotationMatrixY(double angle) {
        // Create rotation using the exponential map with axis (0,1,0)
        Eigen::Vector3d axis = Eigen::Vector3d::UnitY();
        return SO3d::exp(angle * axis).matrix();
    }

    /*!
     * @brief Computes the rotation matrix around the Z-axis using SO3
     *
     * Uses the Lie group SO3 implementation for better numerical stability
     * and consistency with other transformations.
     *
     * @param angle The rotation angle in radians
     * @return RotMat A 3x3 rotation matrix representing the rotation around the Z-axis
     */
    RotMat rotationMatrixZ(double angle) {
        // Create rotation using the exponential map with axis (0,0,1)
        Eigen::Vector3d axis = Eigen::Vector3d::UnitZ();
        return SO3d::exp(angle * axis).matrix();
    }
    sofa::Data<bool> d_debug;

    using Inherit1::fromModels1;
    using Inherit1::fromModels2;
    using Inherit1::toModels;

    sofa::core::State<In1>*m_strain_state;
    sofa::core::State<In2>*m_rigid_base;
    sofa::core::State<Out>*m_global_frames;

protected:
    /// Constructor
    BaseCosseratMapping();

    /// Destructor
    ~BaseCosseratMapping() override = default;

    /**
     * @brief Computes the exponential map for SE(3) using Lie group theory
     * 
     * This function calculates the frame transformation resulting from applying
     * the exponential map to a twist vector scaled by the section length.
     * 
     * @param sub_section_length The length of the beam section
     * @param k The twist vector (angular and linear velocity)
     * @param frame_i The resulting frame transformation
     */
    void computeExponentialSE3(const double &sub_section_length,
                               const Coord1 &k, Frame &frame_i);

    /**
     * @brief Computes the adjoint matrix of a transformation frame
     * 
     * The adjoint matrix is used to transform twists between different reference frames.
     * 
     * @param frame The transformation frame
     * @param adjoint Output adjoint matrix
     */
    void computeAdjoint(const Frame &frame, TangentTransform &adjoint);
    
    /**
     * @brief Computes the adjoint matrix from a 6D vector representation
     * 
     * @param twist The twist vector (angular and linear velocity)
     * @param adjoint Output adjoint matrix
     */
    void computeAdjoint(const Vec6 &twist, Mat6x6 &adjoint);

    void computeCoAdjoint(const Frame &frame, Mat6x6 &coAdjoint);

    void updateExponentialSE3(const vector<Coord1> &strain_state);
    void updateTangExpSE3(const vector<Coord1> &inDeform);

    void computeTangExp(double &x, const Coord1 &k, Mat6x6 &TgX);
    void computeTangExpImplementation(double &x, const Vec6 &k, Mat6x6 &TgX);

    [[maybe_unused]] Vec6
    computeETA(const Vec6 &baseEta, const vector<Deriv1> &k_dot, double abs_input);
    Mat4x4 computeLogarithm(const double &x, const Mat4x4 &gX);
};

#if !defined(SOFA_COSSERAT_CPP_BaseCosseratMapping)
extern template class SOFA_COSSERAT_API
BaseCosseratMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
                 sofa::defaulttype::Rigid3Types>;
extern template class SOFA_COSSERAT_API
BaseCosseratMapping<sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types,
                 sofa::defaulttype::Rigid3Types>;
#endif

} // namespace cosserat::mapping
