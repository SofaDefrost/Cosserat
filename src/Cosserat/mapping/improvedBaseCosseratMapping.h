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

// Liegroups headers
#include <Cosserat/liegroups/Types.h>
#include <Cosserat/liegroups/SO3.h>
#include <Cosserat/liegroups/SE3.h>
#include <Cosserat/liegroups/Utils.h>

#include <sofa/core/Multi2Mapping.h>

namespace Cosserat::mapping
{

// Common type aliases for cleaner code
using namespace sofa::type;
using namespace sofa::component::cosserat::liegroups;

// Type aliases for Lie groups and common types
using TypesD = Types<double>;
using SO3d = sofa::component::cosserat::liegroups::SO3<double>;
using SE3d = sofa::component::cosserat::liegroups::SE3<double>;

// Legacy types (for backward compatibility)
using SE3Legacy = sofa::type::Matrix4;
using se3Legacy = sofa::type::Matrix4;

// Eigen aliases for cleaner code
using Matrix3d = Eigen::Matrix3d;
using Matrix4d = Eigen::Matrix4d;
using Vector3d = Eigen::Vector3d;
using Vector6d = Eigen::Matrix<double, 6, 1>;

// SOFA type aliases
using Frame = Cosserat::type::Frame;
using TangentTransform = Cosserat::type::TangentTransform;
using RotMat = Cosserat::type::RotMat;

/*!
 * \class BaseCosseratMapping
 * @brief Base class for Cosserat rod mappings in SOFA framework
 *
 * This class provides the foundation for implementing Cosserat rod mappings,
 * which are used to map between different representations of a Cosserat rod's
 * configuration and deformation.
 *
 * The implementation uses modern Lie group theory through the LieGroupBase
 * hierarchy to provide numerically stable operations on SE(3) and SO(3).
 *
 * @tparam TIn1 The first input type for the mapping (typically strain state)
 * @tparam TIn2 The second input type for the mapping (typically rigid base)
 * @tparam TOut The output type for the mapping (typically rigid frames)
 */
template <class TIn1, class TIn2, class TOut>
class BaseCosseratMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut>
{
public:
    SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE3(BaseCosseratMapping, TIn1, TIn2, TOut),
                        SOFA_TEMPLATE3(sofa::core::Multi2Mapping,TIn1, TIn2, TOut));

    // Input and output typedefs
    typedef TIn1 In1;
    typedef TIn2 In2;
    typedef TOut Out;

    // Coordinates and derivatives
    using Coord1 = sofa::Coord_t<In1>;
    using Deriv1 = sofa::Deriv_t<In1>;
    using OutCoord = sofa::Coord_t<Out>;

    // Matrix/vector typedefs for readability
    using Matrix6d = TypesD::Matrix<6, 6>;
    using Matrix4d = TypesD::Matrix<4, 4>;
    using Matrix3d = TypesD::Matrix<3, 3>;
    using Vector6d = TypesD::Vector<6>;
    using Vector3d = TypesD::Vector<3>;

    /*===========COSSERAT VECTORS ======================*/
    unsigned int m_indexInput;
    vector<OutCoord> m_vecTransform;

    vector<Frame> m_framesExponentialSE3Vectors;
    vector<Frame> m_nodesExponentialSE3Vectors;
    vector<Matrix4d> m_nodesLogarithmSE3Vectors;

    vector<unsigned int> m_indicesVectors;
    vector<unsigned int> m_indicesVectorsDraw;

    vector<double> m_beamLengthVectors;
    vector<double> m_framesLengthVectors;

    vector<Vector6d> m_nodesVelocityVectors;
    vector<Matrix6d> m_nodesTangExpVectors;
    vector<Matrix6d> m_framesTangExpVectors;
    vector<Vector6d> m_totalBeamForceVectors;

    vector<Matrix6d> m_nodeAdjointVectors;

    /**
     * @brief Compute the co-adjoint matrix of a transformation frame
     * 
     * The co-adjoint matrix is the transpose of the adjoint matrix, used
     * to transform wrenches (force-torque pairs) between coordinate frames.
     * 
     * @param frame The transformation frame
     * @param coAdjoint Output co-adjoint matrix
     */
    void computeCoAdjoint(const Frame& frame, Matrix6d& coAdjoint) {
        // Convert Frame to SE3d
        SE3d se3 = frameToSE3(frame);
        
        // Get the adjoint matrix and transpose it
        Matrix6d adjoint = se3.adjoint();
        coAdjoint = adjoint.transpose();
    }

    /**
     * @brief Update exponential vectors for all frames and nodes
     * 
     * @param strain_state Vector of strain states
     */
    void updateExponentialSE3(const vector<Coord1>& strain_state);
    
    /**
     * @brief Update tangent exponential vectors
     * 
     * @param inDeform Vector of deformations
     */
    void updateTangExpSE3(const vector<Coord1>& inDeform);

    /**
     * @brief Compute tangent exponential map
     * 
     * This function computes the right-trivialized tangent of the exponential map,
     * which is useful for calculating Jacobians in Lie group settings.
     * 
     * @param x Parameter for tangent map
     * @param k Strain vector
     * @param TgX Output tangent matrix
     */
    void computeTangExp(double& x, const Coord1& k, Matrix6d& TgX) {
        // Convert to Vector6d if needed then call implementation
        Vector6d kVec;
        for(int i = 0; i < 6; ++i) {
            kVec(i) = k[i];
        }
        computeTangExpImplementation(x, kVec, TgX);
    }
    
    /**
     * @brief Implementation of tangent exponential map
     * 
     * Uses the SE3 dexp function to compute the right-trivialized tangent
     * of the exponential map.
     * 
     * @param x Parameter for tangent map (scaling factor)
     * @param k Strain vector
     * @param TgX Output tangent matrix
     */
    void computeTangExpImplementation(double& x, const Vector6d& k, Matrix6d& TgX) {
        // Use scaled version of the twist vector
        Vector6d scaledK = x * k;
        
        // Use SE3's dexp function (right-trivialized derivative of exp)
        TgX = SE3d::dexp(scaledK);
    }

    /**
     * @brief Compute eta vector for a given input
     * 
     * Integrates the strain rate along the rod to get the body velocity.
     * 
     * @param baseEta Base eta vector
     * @param k_dot Vector of strain rates
     * @param abs_input Position along the rod
     * @return Vector6d Computed eta vector
     */
    [[maybe_unused]] Vector6d
    computeETA(const Vector6d& baseEta, const vector<Deriv1>& k_dot, double abs_input);
    
    /**
     * @brief Compute logarithm map using SE3
     * 
     * Converts a transformation matrix to its corresponding Lie algebra element.
     * 
     * @param x Scaling factor
     * @param gX Transformation matrix
     * @return Matrix4d Logarithm of the transformation
     */
    Matrix4d computeLogarithm(const double& x, const Matrix4d& gX) {
        // Convert Matrix4d to SE3d
        Eigen::Matrix4d eigenMat;
        for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < 4; ++j) {
                eigenMat(i, j) = gX[i][j];
            }
        }
        
        SE3d se3;
        se3.fromMatrix(eigenMat);
        
        // Compute log and scale
        Vector6d logVec = se3.log();
        logVec *= x;
        
        // Convert to matrix form using hat operator
        return SE3d::hat(logVec);
    }

    /**
     * @brief Computes the adjoint matrix from a twist vector
     * 
     * Creates the adjoint representation for the Lie algebra element.
     * 
     * @param twist The twist vector (angular and linear velocity)
     * @param adjoint Output adjoint matrix
     */
    void computeAdjoint(const Vector6d& twist, Matrix6d& adjoint) {
        // Use SE3's ad function (adjoint of the Lie algebra element)
        adjoint = SE3d::ad(twist);
    }
    
    /**
     * @brief Updates velocity state using Lie group operations
     * 
     * Implements proper velocity updates using adjoint transformations and
     * the Lie group functionality to propagate velocities along the beam.
     * 
     * @param k_dot Strain rates (angular and linear velocity derivatives)
     * @param base_velocity Base node velocity in body coordinates
     */
    void updateVelocityState(const vector<Deriv1>& k_dot, const Vector6d& base_velocity);
    
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
        const Vector6d& source_velocity,
        const Frame& target_frame,
        Vector6d& target_velocity) {
        
        // Convert frames to SE3
        SE3d source_se3 = frameToSE3(source_frame);
        SE3d target_se3 = frameToSE3(target_frame);
        
        // Compute the relative transformation
        SE3d rel_transform = target_se3.inverse() * source_se3;
        
        // Transform velocity using adjoint
        target_velocity = rel_transform.adjoint() * source_velocity;
    }

    /**
     * @brief Compute the angle parameter for logarithm calculation
     *
     * Extracts the rotation angle from a transformation matrix.
     *
     * @param x Scaling factor
     * @param gX Transformation matrix
     * @return double The angle parameter
     */
    double computeTheta(const double& x, const Matrix4d& gX) {
        // Extract rotation part
        Matrix3d R;
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                R(i, j) = gX[i][j];
            }
        }
        
        // Convert to SO3 and get the angle from logarithm
        SO3d so3;
        so3.fromMatrix(R);
        Vector3d logVec = so3.log();
        
        // Return scaled angle
        return x * logVec.norm();
    }

    /**
     * @brief Extract rotation matrix from a Frame using SO3
     * 
     * @param frame The input Frame containing orientation as a quaternion
     * @return Matrix3d The 3x3 rotation matrix
     */
    Matrix3d extractRotMatrix(const Frame& frame) {
        // Convert quaternion to SO3
        SO3d so3;
        so3.fromQuaternion(Eigen::Quaterniond(frame.getOrientation()[3],
                                             frame.getOrientation()[0],
                                             frame.getOrientation()[1],
                                             frame.getOrientation()[2]));
        return so3.matrix();
    }

    /**
     * @brief Build a projector matrix from a Frame
     *
     * @param T The transformation frame
     * @return TangentTransform The projector matrix
     */
    TangentTransform buildProjector(const Frame& T) {
        // Convert to SE3
        SE3d se3 = frameToSE3(T);
        
        // The projector is essentially the adjoint matrix
        Matrix6d adjoint = se3.adjoint();
        
        // Convert to TangentTransform format
        TangentTransform projector;
        for(int i = 0; i < 6; ++i) {
            for(int j = 0; j < 6; ++j) {
                projector[i][j] = adjoint(i, j);
            }
        }
        
        return projector;
    }

    /**
     * @brief Create a skew-symmetric matrix from a vector using SO3::hat
     * 
     * @param u The input 3D vector
     * @return sofa::type::Matrix3 The skew-symmetric matrix
     */
    sofa::type::Matrix3 getTildeMatrix(const Vec3& u) {
        // Use SO3's hat operator
        Eigen::Vector3d v(u[0], u[1], u[2]);
        Eigen::Matrix3d skew = SO3d::hat(v);
        
        // Convert to SOFA Matrix3
        sofa::type::Matrix3 result;
        for(int i = 0; i < 3; ++i) {
            for(int j = 0; j < 3; ++j) {
                result[i][j] = skew(i, j);
            }
        }
        
        return result;
    }

    /**
     * @brief Print a matrix using modern logging
     */
    void printMatrix(const Matrix6d& matrix);
    
    /**
     * @brief Convert a SOFA Frame to an SE3 transformation
     * 
     * @param frame The input SOFA Frame
     * @return SE3d The SE3 transformation
     */
    SE3d frameToSE3(const Frame& frame) {
        // Extract quaternion and position
        Eigen::Quaterniond quat(frame.getOrientation()[3],
                               frame.getOrientation()[0],
                               frame.getOrientation()[1],
                               frame.getOrientation()[2]);
        
        Eigen::Vector3d pos(frame.getCenter()[0],
                           frame.getCenter()[1],
                           frame.getCenter()[2]);
                           
        // Create SE3 from rotation and translation
        SO3d rotation;
        rotation.fromQuaternion(quat);
        
        return SE3d(rotation, pos);
    }
    
    /**
     * @brief Convert an SE3 transformation to a SOFA Frame
     * 
     * @param transform The input SE3 transformation
     * @return Frame The SOFA Frame
     */
    Frame SE3ToFrame(const SE3d& transform) {
        // Extract rotation as quaternion
        SO3d rot = transform.rotation();
        Eigen::Quaterniond quat = rot.toQuaternion();
        
        // Extract translation
        Eigen::Vector3d trans = transform.translation();
        
        // Create Frame
        Frame result;
        result.getOrientation()[0] = quat.x();
        result.getOrientation()[1] = quat.y();
        result.getOrientation()[2] = quat.z();
        result.getOrientation()[3] = quat.w();
        
        result.getCenter()[0] = trans(0);
        result.getCenter()[1] = trans(1);
        result.getCenter()[2] = trans(2);
        
        return result;
    }

    Matrix4d convertTransformToMatrix4x4(const Frame& T) {
        SE3d se3 = frameToSE3(T);
        return se3.matrix();
    }
    
    Vector6d piecewiseLogmap(const Matrix4d& g_x) {
        // Convert to SE3 and compute logarithm
        SE3d se3;
        se3.fromMatrix(g_x);
        return se3.log();
    }

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
        Vector3d axis = Vector3d::UnitX();
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
        Vector3d axis = Vector3d::UnitY();
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
        Vector3d axis = Vector3d::UnitZ();
        return SO3d::exp(angle * axis).matrix();
    }
    
    sofa::Data<bool> d_debug;

    using Inherit1 = sofa::core::Multi2Mapping<TIn1, TIn2, TOut>;
    using Inherit1::fromModels1;
    using Inherit1::fromModels2;
    using Inherit1::toModels;

    sofa::core::State<In1>* m_strain_state;
    sofa::core::State<In2>* m_rigid_base;
    sofa::core::State<Out>* m_global_frames;

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
    void computeExponentialSE3(const double& sub_section_length,
                               const Coord1& k, Frame& frame_i) {
        // Convert Coord1 to Vector6d
        Vector6d twist;
        for(int i = 0; i < 6; ++i) {
            twist(i) = k[i];
        }
        
        // Scale by section length
        twist *= sub_section_length;
        
        // Compute exponential map
        SE3d transform = SE3d::exp(twist);
        
        // Convert to Frame
        frame_i = SE3ToFrame(transform);
    }

    /**
     * @brief Computes the adjoint matrix of a transformation frame
     * 
     * The adjoint matrix is used to transform twists between different reference frames.
     * 
     * @param frame The transformation frame
     * @param adjoint Output adjoint matrix
     */
    void computeAdjoint(const Frame& frame, TangentTransform& adjoint) {
        // Convert Frame to SE3
        SE3d se3 = frameToSE3(frame);
        
        // Get the adjoint matrix
        Matrix6d adjMat = se3.adjoint();
        
        // Convert to TangentTransform
        for(int i = 0; i < 6; ++i) {
            for(int j = 0; j < 6; ++j) {
                adjoint[i][j] = adjMat(i, j);
            }
        }
    }
};

#if !defined(SOFA_COSSERAT_CPP_BaseCosseratMapping)
extern template class SOFA_COSSERAT_API
BaseCosseratMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
                 sofa::defaulttype::Rigid3Types>;
extern template class SOFA_COSSERAT_API
BaseCosseratMapping<sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types,
                 sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
