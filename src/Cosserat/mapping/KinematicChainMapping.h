/******************************************************************************
 * Cosserat plugin — KinematicChainMapping
 *
 * Maps joint coordinates q (Vec1d, one per revolute/prismatic joint) to the world poses of the segments of a
 * kinematic tree (Rigid3d), by product of exponentials on SE(3) (src/liegroups/KinematicChain.h).
 *
 * It plays the role of SOFA's ArticulatedSystemMapping, with a different description of the joints:
 *   - one entry per segment in jointTypes / axes / pivots / restPoses / parents instead of ArticulationCenter and
 *     Articulation components; the rest pose is a full pose (position AND orientation), the axis and the pivot are
 *     given in the frame of the parent segment, so there is a single unambiguous convention;
 *   - a segment hangs from the base or from any earlier segment: trees are supported, serial chains are the default;
 *   - input2 (optional) is the pose of a moving base (Rigid3d), as in ArticulatedSystemMapping.
 *
 * Velocities (applyJ) and forces (applyJT) use SOFA's Rigid3d convention: linear velocity of the segment origin and
 * angular velocity, both in the world frame; forces and torques (about the segment origin) in the world frame.
 *
 * Dynamics support (what a mass defined on the output frames needs to give exact joint-space dynamics):
 *   - applyDJT: geometric stiffness d(J^T)/dx . dx . f for the constant child force f (data geometricStiffness);
 *   - computeBiasAcceleration / computeAccFromMapping: the velocity-dependent part of the frame accelerations, d/dt(J dx)
 *     for zero generalised accelerations. SOFA does not apply it by itself: KinematicChainBiasForceField turns it into the
 *     inertial force -M (Jdot qdot) on the frames, which completes  (J^T M J) qddot = J^T (f - M Jdot qdot).
 * The maths lives in liegroups/KinematicChain.h (biasAcceleration, geometricStiffness) and is tested against finite
 * differences in simulations/articulated/lie_tests/test_lie_chain_derivatives.py.
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>

#include <sofa/core/Multi2Mapping.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Quat.h>
#include <sofa/type/vector.h>

#include <liegroups/KinematicChain.h>

#include <string>
#include <vector>

namespace Cosserat::mapping {

template<class TIn1, class TIn2, class TOut>
class KinematicChainMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut> {
public:
    SOFA_CLASS(SOFA_TEMPLATE3(KinematicChainMapping, TIn1, TIn2, TOut),
               SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut));
    typedef sofa::core::Multi2Mapping<TIn1, TIn2, TOut> Inherit;

    typedef TIn1 In1;
    typedef TIn2 In2;
    typedef TOut Out;

    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;
    typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef typename In1::Deriv In1Deriv;
    typedef typename In1::Real In1Real;
    typedef sofa::Data<In1VecCoord> In1DataVecCoord;
    typedef sofa::Data<In1VecDeriv> In1DataVecDeriv;
    typedef sofa::Data<In1MatrixDeriv> In1DataMatrixDeriv;

    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;
    typedef typename In2::MatrixDeriv In2MatrixDeriv;
    typedef typename In2::Coord In2Coord;
    typedef typename In2::Deriv In2Deriv;
    typedef sofa::Data<In2VecCoord> In2DataVecCoord;
    typedef sofa::Data<In2VecDeriv> In2DataVecDeriv;
    typedef sofa::Data<In2MatrixDeriv> In2DataMatrixDeriv;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef sofa::Data<OutVecCoord> OutDataVecCoord;
    typedef sofa::Data<OutVecDeriv> OutDataVecDeriv;
    typedef sofa::Data<OutMatrixDeriv> OutDataMatrixDeriv;

    using Chain = sofa::component::cosserat::liegroups::KinematicChain<double>;
    using SE3d = sofa::component::cosserat::liegroups::SE3<double>;
    using Vec3 = sofa::type::Vec3d;

    /// One entry per segment (shorter lists are padded with the defaults noted below).
    sofa::Data<sofa::type::vector<std::string>> d_jointTypes; ///< "revolute", "prismatic" or "fixed"
    sofa::Data<sofa::type::vector<Vec3>> d_axes;              ///< joint axis in the parent frame (default 0 0 1)
    sofa::Data<sofa::type::vector<Vec3>> d_pivots;            ///< a point on the revolute axis, parent frame (default 0)
    sofa::Data<sofa::type::vector<OutCoord>> d_restPoses;     ///< pose of the segment in its parent at q = 0 (default identity)
    sofa::Data<sofa::type::vector<int>> d_parents;            ///< parent segment, -1 = base; empty = serial chain
    sofa::Data<unsigned int> d_baseIndex;                     ///< index of the base pose in input2
    sofa::Data<bool> d_geometricStiffness;                    ///< applyDJT: add the geometric stiffness (default true)
    sofa::Data<bool> d_drawJoints;                            ///< draw the joint axes (needs showMappings)
    sofa::Data<double> d_drawSize;                            ///< length of the drawn axes

protected:
    KinematicChainMapping();
    ~KinematicChainMapping() override = default;

public:
    void init() override;
    void reinit() override;
    void draw(const sofa::core::visual::VisualParams *vparams) override;

    using Inherit::apply;
    using Inherit::applyJ;
    using Inherit::applyJT;

    void apply(const sofa::core::MechanicalParams *, const sofa::type::vector<OutDataVecCoord *> &dataVecOutPos,
               const sofa::type::vector<const In1DataVecCoord *> &dataVecIn1Pos,
               const sofa::type::vector<const In2DataVecCoord *> &dataVecIn2Pos) override;

    void applyJ(const sofa::core::MechanicalParams *, const sofa::type::vector<OutDataVecDeriv *> &dataVecOutVel,
                const sofa::type::vector<const In1DataVecDeriv *> &dataVecIn1Vel,
                const sofa::type::vector<const In2DataVecDeriv *> &dataVecIn2Vel) override;

    void applyJT(const sofa::core::MechanicalParams *, const sofa::type::vector<In1DataVecDeriv *> &dataVecOut1Force,
                 const sofa::type::vector<In2DataVecDeriv *> &dataVecOut2Force,
                 const sofa::type::vector<const OutDataVecDeriv *> &dataVecInForce) override;

    void applyJT(const sofa::core::ConstraintParams *, const sofa::type::vector<In1DataMatrixDeriv *> &dataMatOut1,
                 const sofa::type::vector<In2DataMatrixDeriv *> &dataMatOut2,
                 const sofa::type::vector<const OutDataMatrixDeriv *> &dataMatIn) override;

    /// Geometric stiffness: parentForce += kFactor * d(J^T)/dx . dx . childForce  (dx = displacement read from the inputs).
    void applyDJT(const sofa::core::MechanicalParams *mparams, sofa::core::MultiVecDerivId parentForce,
                  sofa::core::ConstMultiVecDerivId childForce) override;

    /// Frame accelerations for the input accelerations: J qddot + (velocity-dependent term).
    void computeAccFromMapping(const sofa::core::MechanicalParams *, const sofa::type::vector<OutDataVecDeriv *> &dataVecOutAcc,
                               const sofa::type::vector<const In1DataVecDeriv *> &dataVecIn1Vel,
                               const sofa::type::vector<const In2DataVecDeriv *> &dataVecIn2Vel,
                               const sofa::type::vector<const In1DataVecDeriv *> &dataVecIn1Acc,
                               const sofa::type::vector<const In2DataVecDeriv *> &dataVecIn2Acc) override;
    using Inherit::computeAccFromMapping;

    /// Velocity-dependent part of the frame accelerations (angular, linear of the origin; world frame), computed from the
    /// current positions (last apply) and the velocities of the input states. Zero accelerations of the inputs assumed.
    void computeBiasAcceleration(OutVecDeriv &bias);

    const sofa::linearalgebra::BaseMatrix *getJ() override { return nullptr; }

    /// The kinematic tree built from the Data (read-only access, e.g. for tests).
    const Chain &chain() const { return m_chain; }

private:
    void buildChain();

    static SE3d toSE3(const OutCoord &c);
    static OutCoord fromSE3(const SE3d &T);

    /// Fills the cached world poses / screws from q and the optional base pose.
    void update(const In1VecCoord &q, const In2VecCoord *base);

    /// J applied to (dq, dbase), shared by applyJ and computeAccFromMapping.
    void mapVelocity(const In1VecDeriv &dq, const In2VecDeriv *dbase, OutVecDeriv &out) const;
    /// Velocity-dependent acceleration for the velocities (dq, dbase) (zero accelerations), added to @p out.
    void addBias(const In1VecDeriv &dq, const In2VecDeriv *dbase, OutVecDeriv &out) const;

    Chain m_chain;
    bool m_valid = false;

    // cache refreshed by apply() and used by applyJ / applyJT
    std::vector<Vec3> m_position;                                 ///< world position of every segment origin
    std::vector<Eigen::Vector3d> m_positionE;                     ///< same, Eigen (for KinematicChain)
    std::vector<SE3d> m_poses;                                    ///< world pose of every segment
    Eigen::Matrix<double, 6, Eigen::Dynamic> m_screws;            ///< world-frame unit screws [w ; r] of every joint
    Vec3 m_baseCenter = Vec3(0, 0, 0);
    bool m_hasBase = false;
};

#if !defined(SOFA_COSSERAT_CPP_KinematicChainMapping)
extern template class SOFA_COSSERAT_API KinematicChainMapping<sofa::defaulttype::Vec1Types, sofa::defaulttype::Rigid3Types,
                                                              sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
