#pragma once
#include <Cosserat/mapping/KinematicChainMapping.h>

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Messaging.h>

#include <optional>

namespace Cosserat::mapping {

using sofa::core::objectmodel::ComponentState;

template<class TIn1, class TIn2, class TOut>
KinematicChainMapping<TIn1, TIn2, TOut>::KinematicChainMapping()
    : d_jointTypes(initData(&d_jointTypes, "jointTypes", "Per segment: revolute, prismatic or fixed")),
      d_axes(initData(&d_axes, "axes", "Per segment: joint axis in the parent frame (default 0 0 1)")),
      d_pivots(initData(&d_pivots, "pivots", "Per segment: a point on the revolute axis, in the parent frame (default 0 0 0)")),
      d_restPoses(initData(&d_restPoses, "restPoses",
                           "Per segment: pose (x y z qx qy qz qw) of the segment in its parent frame at q = 0 (default identity)")),
      d_parents(initData(&d_parents, "parents", "Per segment: index of the parent segment, -1 for the base. Empty = serial chain")),
      d_baseIndex(initData(&d_baseIndex, 0u, "baseIndex", "Index of the base pose in input2 (used only when input2 is given)")),
      d_drawJoints(initData(&d_drawJoints, true, "drawJoints", "Draw the joint axes when showMappings is on")),
      d_drawSize(initData(&d_drawSize, 1.0, "drawSize", "Length of the drawn joint axes")) {}

template<class TIn1, class TIn2, class TOut>
typename KinematicChainMapping<TIn1, TIn2, TOut>::SE3d KinematicChainMapping<TIn1, TIn2, TOut>::toSE3(const OutCoord &c) {
    const auto &o = c.getOrientation();
    Eigen::Quaterniond q(o[3], o[0], o[1], o[2]); // Eigen ctor order: w, x, y, z
    q.normalize();
    return SE3d(q, Eigen::Vector3d(c.getCenter()[0], c.getCenter()[1], c.getCenter()[2]));
}

template<class TIn1, class TIn2, class TOut>
typename KinematicChainMapping<TIn1, TIn2, TOut>::OutCoord KinematicChainMapping<TIn1, TIn2, TOut>::fromSE3(const SE3d &T) {
    const Eigen::Quaterniond q = T.rotation().quaternion();
    const Eigen::Vector3d &t = T.translation();
    OutCoord c;
    c.getCenter() = Vec3(t.x(), t.y(), t.z());
    c.getOrientation() = sofa::type::Quat<double>(q.x(), q.y(), q.z(), q.w());
    return c;
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::buildChain() {
    m_chain = Chain();
    m_valid = false;
    const auto &types = d_jointTypes.getValue();
    const auto &axes = d_axes.getValue();
    const auto &pivots = d_pivots.getValue();
    const auto &rests = d_restPoses.getValue();
    const auto &parents = d_parents.getValue();
    if (types.empty()) {
        msg_error() << "jointTypes is empty: nothing to map.";
        return;
    }
    try {
        for (std::size_t i = 0; i < types.size(); ++i) {
            const Vec3 a = i < axes.size() ? axes[i] : Vec3(0, 0, 1);
            const Vec3 p = i < pivots.size() ? pivots[i] : Vec3(0, 0, 0);
            const SE3d rest = i < rests.size() ? toSE3(rests[i]) : SE3d();
            const int parent = i < parents.size() ? parents[i] : Chain::Previous;
            if (types[i] == "revolute")
                m_chain.addRevolute(Eigen::Vector3d(a[0], a[1], a[2]), Eigen::Vector3d(p[0], p[1], p[2]), rest, parent);
            else if (types[i] == "prismatic")
                m_chain.addPrismatic(Eigen::Vector3d(a[0], a[1], a[2]), rest, parent);
            else if (types[i] == "fixed")
                m_chain.addFixed(rest, parent);
            else {
                msg_error() << "jointTypes[" << i << "] = '" << types[i] << "': expected revolute, prismatic or fixed.";
                return;
            }
        }
    } catch (const std::exception &e) {
        msg_error() << "Invalid kinematic tree: " << e.what();
        return;
    }
    m_valid = true;
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::init() {
    buildChain();
    if (this->getFromModels1().empty() || this->getToModels().empty()) {
        msg_error() << "input1 (joint coordinates, Vec1d) and output (segment poses, Rigid3d) are required.";
        m_valid = false;
    } else {
        if (static_cast<int>(this->getFromModels1()[0]->getSize()) != m_chain.nDof()) {
            msg_error() << "input1 has " << this->getFromModels1()[0]->getSize() << " coordinates but the tree has "
                        << m_chain.nDof() << " joints.";
            m_valid = false;
        }
        if (static_cast<int>(this->getToModels()[0]->getSize()) != m_chain.nSegments()) {
            msg_error() << "output has " << this->getToModels()[0]->getSize() << " frames but the tree has "
                        << m_chain.nSegments() << " segments.";
            m_valid = false;
        }
    }
    this->d_componentState.setValue(m_valid ? ComponentState::Valid : ComponentState::Invalid);
    Inherit::init();
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::reinit() {
    buildChain();
    this->d_componentState.setValue(m_valid ? ComponentState::Valid : ComponentState::Invalid);
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::update(const In1VecCoord &q, const In2VecCoord *base) {
    Chain::VectorX qv(m_chain.nDof());
    for (int k = 0; k < m_chain.nDof(); ++k) qv[k] = q[k][0];

    m_hasBase = base != nullptr && d_baseIndex.getValue() < base->size();
    SE3d T0; // base pose, identity when input2 is absent
    if (m_hasBase) {
        T0 = toSE3((*base)[d_baseIndex.getValue()]);
        m_baseCenter = Vec3(T0.translation().x(), T0.translation().y(), T0.translation().z());
    } else
        m_baseCenter = Vec3(0, 0, 0);

    const std::vector<SE3d> local = m_chain.allPoses(qv);
    std::vector<SE3d> world(local.size());
    m_position.resize(local.size());
    for (std::size_t i = 0; i < local.size(); ++i) {
        world[i] = T0 * local[i];
        m_position[i] = Vec3(world[i].translation().x(), world[i].translation().y(), world[i].translation().z());
    }
    // unit screws seen from the base frame, then carried to the world frame by the base pose
    m_screws = T0.computeAdjoint() * m_chain.jointScrews(local);
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams *,
                                                    const sofa::type::vector<OutDataVecCoord *> &dataVecOutPos,
                                                    const sofa::type::vector<const In1DataVecCoord *> &dataVecIn1Pos,
                                                    const sofa::type::vector<const In2DataVecCoord *> &dataVecIn2Pos) {
    if (!m_valid || dataVecOutPos.empty() || dataVecIn1Pos.empty()) return;
    const In1VecCoord &q = dataVecIn1Pos[0]->getValue();
    const In2VecCoord *base = dataVecIn2Pos.empty() || !dataVecIn2Pos[0] ? nullptr : &dataVecIn2Pos[0]->getValue();
    if (static_cast<int>(q.size()) != m_chain.nDof()) return;
    update(q, base);

    Chain::VectorX qv(m_chain.nDof());
    for (int k = 0; k < m_chain.nDof(); ++k) qv[k] = q[k][0];
    const std::vector<SE3d> local = m_chain.allPoses(qv);
    SE3d T0;
    if (m_hasBase) T0 = toSE3((*base)[d_baseIndex.getValue()]);

    OutVecCoord &out = *dataVecOutPos[0]->beginEdit();
    out.resize(local.size());
    for (std::size_t i = 0; i < local.size(); ++i) out[i] = fromSE3(T0 * local[i]);
    dataVecOutPos[0]->endEdit();
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams *,
                                                     const sofa::type::vector<OutDataVecDeriv *> &dataVecOutVel,
                                                     const sofa::type::vector<const In1DataVecDeriv *> &dataVecIn1Vel,
                                                     const sofa::type::vector<const In2DataVecDeriv *> &dataVecIn2Vel) {
    if (!m_valid || dataVecOutVel.empty() || dataVecIn1Vel.empty()) return;
    const In1VecDeriv &dq = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv *dbase = dataVecIn2Vel.empty() || !dataVecIn2Vel[0] ? nullptr : &dataVecIn2Vel[0]->getValue();

    Vec3 vb(0, 0, 0), wb(0, 0, 0);
    if (m_hasBase && dbase && d_baseIndex.getValue() < dbase->size()) {
        vb = getVCenter((*dbase)[d_baseIndex.getValue()]);
        wb = getVOrientation((*dbase)[d_baseIndex.getValue()]);
    }

    OutVecDeriv &out = *dataVecOutVel[0]->beginEdit();
    out.resize(m_chain.nSegments());
    for (int i = 0; i < m_chain.nSegments(); ++i) {
        Vec3 w = wb;
        Vec3 v = vb + sofa::type::cross(wb, m_position[i] - m_baseCenter);
        for (int k : m_chain.ancestorDofs(i)) {
            const Vec3 wk(m_screws(0, k), m_screws(1, k), m_screws(2, k));
            const Vec3 rk(m_screws(3, k), m_screws(4, k), m_screws(5, k));
            const double rate = dq[k][0];
            w += wk * rate;
            v += (rk + sofa::type::cross(wk, m_position[i])) * rate;
        }
        getVCenter(out[i]) = v;
        getVOrientation(out[i]) = w;
    }
    dataVecOutVel[0]->endEdit();
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::MechanicalParams *,
                                                      const sofa::type::vector<In1DataVecDeriv *> &dataVecOut1Force,
                                                      const sofa::type::vector<In2DataVecDeriv *> &dataVecOut2Force,
                                                      const sofa::type::vector<const OutDataVecDeriv *> &dataVecInForce) {
    if (!m_valid || dataVecOut1Force.empty() || dataVecInForce.empty()) return;
    const OutVecDeriv &f = dataVecInForce[0]->getValue();
    const int n = m_chain.nSegments();

    // wrench of every segment at the world origin, [torque + p x force ; force], accumulated from the leaves up:
    // all the wrenches share the reference point, so a joint sees the plain sum of what hangs below it.
    std::vector<sofa::type::Vec<6, double>> acc(n);
    for (int i = 0; i < n; ++i) {
        const Vec3 force = getVCenter(f[i]);
        const Vec3 torque = getVOrientation(f[i]) + sofa::type::cross(m_position[i], force);
        for (int a = 0; a < 3; ++a) {
            acc[i][a] = torque[a];
            acc[i][3 + a] = force[a];
        }
    }
    sofa::type::Vec<6, double> total;
    for (int i = n - 1; i >= 0; --i) {
        const int p = m_chain.segments()[i].parent;
        if (p >= 0) acc[p] += acc[i];
        else total += acc[i];
    }

    In1VecDeriv &out1 = *dataVecOut1Force[0]->beginEdit();
    for (int i = 0; i < n; ++i) {
        const int k = m_chain.segments()[i].dofIndex;
        if (k < 0) continue;
        double q = 0;
        for (int a = 0; a < 6; ++a) q += m_screws(a, k) * acc[i][a];
        out1[k][0] += static_cast<In1Real>(q);
    }
    dataVecOut1Force[0]->endEdit();

    if (m_hasBase && !dataVecOut2Force.empty() && dataVecOut2Force[0] && d_baseIndex.getValue() < dataVecOut2Force[0]->getValue().size()) {
        In2VecDeriv &out2 = *dataVecOut2Force[0]->beginEdit();
        const Vec3 force(total[3], total[4], total[5]);
        const Vec3 torqueOrigin(total[0], total[1], total[2]);
        getVCenter(out2[d_baseIndex.getValue()]) += force;
        getVOrientation(out2[d_baseIndex.getValue()]) += torqueOrigin - sofa::type::cross(m_baseCenter, force);
        dataVecOut2Force[0]->endEdit();
    }
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::ConstraintParams *,
                                                      const sofa::type::vector<In1DataMatrixDeriv *> &dataMatOut1,
                                                      const sofa::type::vector<In2DataMatrixDeriv *> &dataMatOut2,
                                                      const sofa::type::vector<const OutDataMatrixDeriv *> &dataMatIn) {
    if (!m_valid || dataMatOut1.empty() || dataMatIn.empty()) return;
    In1MatrixDeriv &out1 = *dataMatOut1[0]->beginEdit();
    In2MatrixDeriv *out2 = (m_hasBase && !dataMatOut2.empty() && dataMatOut2[0]) ? dataMatOut2[0]->beginEdit() : nullptr;
    const OutMatrixDeriv &in = dataMatIn[0]->getValue();

    for (auto rowIt = in.begin(); rowIt != in.end(); ++rowIt) {
        auto colIt = rowIt.begin();
        if (colIt == rowIt.end()) continue;
        auto o1 = out1.writeLine(rowIt.index());
        std::optional<typename In2MatrixDeriv::RowIterator> o2;
        if (out2) o2.emplace(out2->writeLine(rowIt.index()));

        for (; colIt != rowIt.end(); ++colIt) {
            const int i = colIt.index();
            const OutDeriv &value = colIt.val();
            const Vec3 v = getVCenter(value);
            const Vec3 torqueOrigin = getVOrientation(value) + sofa::type::cross(m_position[i], v);
            for (int k : m_chain.ancestorDofs(i)) {
                const double coef = m_screws(0, k) * torqueOrigin[0] + m_screws(1, k) * torqueOrigin[1] +
                                    m_screws(2, k) * torqueOrigin[2] + m_screws(3, k) * v[0] + m_screws(4, k) * v[1] +
                                    m_screws(5, k) * v[2];
                In1Deriv d;
                d[0] = static_cast<In1Real>(coef);
                o1.addCol(k, d);
            }
            if (out2) {
                In2Deriv db;
                getVCenter(db) = v;
                getVOrientation(db) = getVOrientation(value) + sofa::type::cross(m_position[i] - m_baseCenter, v);
                o2->addCol(d_baseIndex.getValue(), db);
            }
        }
    }
    dataMatOut1[0]->endEdit();
    if (out2) dataMatOut2[0]->endEdit();
}

template<class TIn1, class TIn2, class TOut>
void KinematicChainMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
    if (!m_valid || !vparams->displayFlags().getShowMappings() || !d_drawJoints.getValue()) return;
    const double s = d_drawSize.getValue();
    std::vector<sofa::type::Vec3> axes, points;
    for (int i = 0; i < m_chain.nSegments(); ++i) {
        const int k = m_chain.segments()[i].dofIndex;
        if (k < 0) continue;
        const Vec3 w(m_screws(0, k), m_screws(1, k), m_screws(2, k));
        const Vec3 r(m_screws(3, k), m_screws(4, k), m_screws(5, k));
        // a point on the joint axis: for a revolute joint the foot of the perpendicular from the origin, w x r
        const bool revolute = m_chain.segments()[i].type == Chain::JointType::Revolute;
        const Vec3 a = revolute ? sofa::type::cross(w, r) : m_position[i];
        const Vec3 dir = revolute ? w : r;
        points.emplace_back(a[0], a[1], a[2]);
        axes.emplace_back(a[0], a[1], a[2]);
        axes.emplace_back(a[0] + s * dir[0], a[1] + s * dir[1], a[2] + s * dir[2]);
    }
    vparams->drawTool()->drawPoints(points, 8, sofa::type::RGBAColor(1, 0.5, 0.5, 1));
    vparams->drawTool()->drawLines(axes, 2, sofa::type::RGBAColor::blue());
}

} // namespace Cosserat::mapping
