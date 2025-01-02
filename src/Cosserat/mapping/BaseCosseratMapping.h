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

using SE3 = sofa::type::Matrix4; ///< The "coordinate" in SE3
using se3 = sofa::type::Matrix4; ///< The "speed" of change of SE3.
using _se3 = Eigen::Matrix4d;
using _SE3 = Eigen::Matrix4d;

using Cosserat::type::Frame;
using Cosserat::type::TangentTransform;
using Cosserat::type::RotMat;


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

    // TODO(dmarchal:2024/06/07): explain why these attributes are unused
    // : yadagolo: Need for the dynamic function, which is not working yet. But the component is in this folder
    // : dmarchal: don't add something that will be used "one day"
    // : dmarchal: it look like as if you should be working in a branch for making new feature and merge it when it is ready.
    [[maybe_unused]] vector<Mat6x6> m_nodeAdjointEtaVectors;
    [[maybe_unused]] vector<Mat6x6> m_frameAdjointEtaVectors;
    [[maybe_unused]] vector<Mat6x6> m_node_coAdjointEtaVectors;
    [[maybe_unused]] vector<Mat6x6> m_frame_coAdjointEtaVectors;

public:
    /********************** Inhertited from BaseObject   **************/
    void init() final override;
    virtual void doBaseCosseratInit() = 0;

    // This function is called by a callback function, which is not the case
    // of the init function
    void update_geometry_info();

    double computeTheta(const double &x, const Mat4x4 &gX);

    static void printMatrix(const Mat6x6 R);

    static sofa::type::Mat3x3 extractRotMatrix(const Frame &frame);

    static TangentTransform buildProjector(const Frame &T);

    static Mat3x3 getTildeMatrix(const Vec3 &u);

    static void buildAdjoint(const Mat3x3 &A, const Mat3x3 &B, Mat6x6 &Adjoint);
    void buildCoAdjoint(const Mat3x3 &A, const Mat3x3 &B, Mat6x6 &coAdjoint);

    static Mat4x4 convertTransformToMatrix4x4(const Frame &T);

    static Vec6 piecewiseLogmap(const _SE3 &g_x);

    /*!
     * @brief Computes the rotation matrix around the X-axis
     *
     * @param angle The rotation angle in radians
     * @return RotMat A 3x3 rotation matrix representing the rotation around the X-axis
     */
    RotMat rotationMatrixX(double angle) {
        Eigen::Matrix3d rotation;
        rotation << 1, 0, 0, 0, cos(angle), -sin(angle), 0, sin(angle), cos(angle);
        return rotation;
    }

    /*!
     * @brief Computes the rotation matrix around the Y-axis
     *
     * @param angle The rotation angle in radians
     * @return RotMat A 3x3 rotation matrix representing the rotation around the Y-axis
     */
    // function... it shouldn't be (re)implemented in a base classe.
    RotMat rotationMatrixY(double angle) {
        Eigen::Matrix3d rotation;
        rotation << cos(angle), 0, sin(angle), 0, 1, 0, -sin(angle), 0, cos(angle);
        return rotation;
    }

    // TODO(dmarchal: 2024/06/07), this looks like a very common utility
    // function... it shouldn't be (re)implemented in a base classe. the type of
    // the data return should also be unified between rotationMatrixX, Y and Z
    static RotMat rotationMatrixZ(const SReal angle) {
        RotMat rotation;
        rotation << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
        return rotation;
    }

protected:
    sofa::Data<vector<double>> d_curv_abs_section;
    sofa::Data<vector<double>> d_curv_abs_frames;
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

    void computeExponentialSE3(const double &sub_section_length,
                               const Coord1 &k,Frame &frame_i);

    // TODO(dmarchal: 2024/06/07):
    //   - clarify the difference between computeAdjoing and buildAdjoint ...
    //   - clarify why we need Transform and Vec6 and TangentTransform & Mat6x6
    void computeAdjoint(const Frame &frame, TangentTransform &adjoint);
    void computeAdjoint(const Vec6 &frame, Mat6x6 &adjoint);

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
