/******************************************************************************
 * Cosserat plugin — RigidBodyMass
 *
 * Full rigid-body mass for Rigid3d states: mass, centre of mass anywhere in the body frame and a
 * general (anisotropic) inertia tensor. SOFA's UniformMass / DiagonalMass apply a CONSTANT inertia
 * matrix to the world-frame angular quantities and add no gyroscopic force, so they are exact only
 * for an isotropic inertia and a centre of mass on the frame origin. This component is the exact
 * Newton-Euler model of a body whose state is (position of the frame origin, orientation):
 *
 *   world quantities   c = R c_b   (centre of mass relative to the frame origin),   I = R I_b R^T
 *   mass matrix        M = [[ m 1      , -m [c]x          ],        acts on (linear acc. of the
 * origin, angular acc.) [ m [c]x   ,  I - m [c]x [c]x ]]        and returns (force, torque about
 * the origin) force              F = [ m g - m w x (w x c) ;  m c x g - w x (I w) - m c x (w x (w x
 * c)) ]
 *
 * (w = angular velocity in the world frame, which is SOFA's convention for the angular part of a
 * Rigid3d velocity.) Used with KinematicChainMapping and KinematicChainBiasForceField it makes
 * route A (SOFA mass on the mapped frames) equal to the spatial-inertia Newton-Euler dynamics of
 * KinematicChainDynamics (route B) for any body inertia.
 *
 * Tangent (addDForce / buildStiffnessMatrix / buildDampingMatrix), exact:
 *   - dF/d(rotation): the weight torque m c x g and the centrifugal / gyroscopic terms depend on R
 * through c and I (d c = dtheta x c,  d I w = dtheta x (I w) - I (dtheta x w));
 *   - dF/dw: the gyroscopic terms. This block is not symmetric (it contains the skew part of the
 * gyroscopic coupling), so it is only given when the Data 'gyroDamping' is true; conjugate gradient
 * may then need a non-symmetric solver.
 *
 * Data (one entry per body, or one entry for all of them):
 *   mass            [kg]
 *   centerOfMass    [m] in the body frame
 *   inertia         (Ixx, Iyy, Izz, Ixy, Ixz, Iyz): the entries of the symmetric inertia matrix
 * about the centre of mass, in the body frame, as in a URDF <inertia> tag (off-diagonal entries are
 * the matrix entries themselves)
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <sofa/core/behavior/BaseLocalForceFieldMatrix.h>
#include <sofa/core/behavior/BaseLocalMassMatrix.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

#include <Eigen/Dense>

namespace Cosserat::mapping {

class SOFA_COSSERAT_API RigidBodyMass
    : public sofa::core::behavior::Mass<sofa::defaulttype::Rigid3Types> {
   public:
    SOFA_CLASS(RigidBodyMass,
               SOFA_TEMPLATE(sofa::core::behavior::Mass, sofa::defaulttype::Rigid3Types));
    using DataTypes = sofa::defaulttype::Rigid3Types;
    using Inherit = sofa::core::behavior::Mass<DataTypes>;
    using Coord = DataTypes::Coord;
    using Deriv = DataTypes::Deriv;
    using VecCoord = DataTypes::VecCoord;
    using VecDeriv = DataTypes::VecDeriv;
    using DataVecCoord = sofa::Data<VecCoord>;
    using DataVecDeriv = sofa::Data<VecDeriv>;
    using Vec3 = Eigen::Vector3d;
    using Mat3 = Eigen::Matrix3d;
    using Vec6 = Eigen::Matrix<double, 6, 1>;
    using Mat6 = Eigen::Matrix<double, 6, 6>;

    sofa::Data<sofa::type::vector<double>> d_mass;
    sofa::Data<sofa::type::vector<sofa::type::Vec3d>> d_centerOfMass;
    sofa::Data<sofa::type::vector<sofa::type::Vec<6, double>>> d_inertia;
    sofa::Data<bool> d_stiffness;    ///< tangent: derivative of the force w.r.t. the orientation
                                     ///< (weight torque, gyroscopic)
    sofa::Data<bool> d_gyroDamping;  ///< tangent: derivative of the gyroscopic force w.r.t. the
                                     ///< angular velocity (not symmetric)

   protected:
    RigidBodyMass()
        : d_mass(initData(&d_mass, sofa::type::vector<double>{1.0}, "mass",
                          "mass of each body [kg] (one value = all bodies)")),
          d_centerOfMass(initData(
              &d_centerOfMass, sofa::type::vector<sofa::type::Vec3d>{sofa::type::Vec3d(0, 0, 0)},
              "centerOfMass", "centre of mass in the body frame [m] (one value = all bodies)")),
          d_inertia(initData(&d_inertia,
                             sofa::type::vector<sofa::type::Vec<6, double>>{
                                 sofa::type::Vec<6, double>(1, 1, 1, 0, 0, 0)},
                             "inertia",
                             "inertia about the centre of mass, body frame: Ixx Iyy Izz Ixy Ixz "
                             "Iyz (matrix entries)")),
          d_stiffness(initData(&d_stiffness, true, "stiffness",
                               "tangent: derivative of the force w.r.t. the orientation")),
          d_gyroDamping(initData(&d_gyroDamping, false, "gyroDamping",
                                 "tangent: derivative of the gyroscopic force w.r.t. the angular "
                                 "velocity (not symmetric)")) {}

   public:
    void init() override {
        Inherit::init();
        if (!this->mstate) {
            msg_error() << "needs a Rigid3d mechanical state in its node.";
            this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
            return;
        }
        reinit();
    }

    void reinit() override {
        const std::size_t n = this->mstate ? this->mstate->getSize() : 0;
        for (std::size_t sz : {d_mass.getValue().size(), d_centerOfMass.getValue().size(),
                               d_inertia.getValue().size()}) {
            if (sz != 1 && sz != n) {
                msg_error() << "mass, centerOfMass and inertia need 1 or " << n
                            << " entries (one per body).";
                this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
                return;
            }
        }
        m_n = n;
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }

    // ---- model
    // -------------------------------------------------------------------------------------------------
    struct Body {
        double m = 1.0;
        Vec3 c = Vec3::Zero();  ///< centre of mass in the world frame, relative to the frame origin
        Mat3 I = Mat3::Identity();  ///< inertia about the centre of mass, world frame
    };

    static Mat3 hat(const Vec3 &a) {
        Mat3 h;
        h << 0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0;
        return h;
    }

    static Mat3 rotation(const Coord &x) {
        const auto q = x.getOrientation();
        return Eigen::Quaterniond(q[3], q[0], q[1], q[2]).normalized().toRotationMatrix();
    }

    Body body(std::size_t i, const Coord &x) const {
        const auto &m = d_mass.getValue();
        const auto &c = d_centerOfMass.getValue();
        const auto &in = d_inertia.getValue();
        const std::size_t k = m.size() == 1 ? 0 : i, kc = c.size() == 1 ? 0 : i,
                          ki = in.size() == 1 ? 0 : i;
        const Mat3 R = rotation(x);
        Mat3 Ib;
        const auto &e = in[ki];
        Ib << e[0], e[3], e[4], e[3], e[1], e[5], e[4], e[5], e[2];
        Body b;
        b.m = m[k];
        b.c = R * Vec3(c[kc][0], c[kc][1], c[kc][2]);
        b.I = R * Ib * R.transpose();
        return b;
    }

    static Mat6 massMatrix(const Body &b) {
        const Mat3 S = hat(b.c);
        Mat6 M;
        M.topLeftCorner<3, 3>() = b.m * Mat3::Identity();
        M.topRightCorner<3, 3>() = -b.m * S;
        M.bottomLeftCorner<3, 3>() = b.m * S;
        M.bottomRightCorner<3, 3>() = b.I - b.m * S * S;
        return M;
    }

    /// Weight and inertial (centrifugal / gyroscopic) wrench about the frame origin, world frame,
    /// for angular velocity w.
    static Vec6 wrench(const Body &b, const Vec3 &w, const Vec3 &g) {
        const Vec3 wwc = w.cross(w.cross(b.c));
        Vec6 F;
        F.head<3>() = b.m * (g - wwc);
        F.tail<3>() = b.m * b.c.cross(g) - w.cross(b.I * w) - b.m * b.c.cross(wwc);
        return F;
    }

    /// d(wrench)/d(rotation) . dtheta   (rotation increment in the world frame: R -> exp(dtheta) R)
    static Vec6 wrenchRotationDerivative(const Body &b, const Vec3 &w, const Vec3 &g,
                                         const Vec3 &dth) {
        const Vec3 dc = dth.cross(b.c);
        const Vec3 Iw = b.I * w;
        const Vec3 dIw = dth.cross(Iw) - b.I * dth.cross(w);  // (d I) w
        const Vec3 wwc = w.cross(w.cross(b.c));
        Vec6 dF;
        dF.head<3>() = -b.m * w.cross(w.cross(dc));
        dF.tail<3>() = b.m * dc.cross(g) - w.cross(dIw) - b.m * dc.cross(wwc) -
                       b.m * b.c.cross(w.cross(w.cross(dc)));
        return dF;
    }

    /// d(wrench)/d(w) . dw
    static Vec6 wrenchVelocityDerivative(const Body &b, const Vec3 &w, const Vec3 &dw) {
        const Vec3 dwwc = dw.cross(w.cross(b.c)) + w.cross(dw.cross(b.c));
        Vec6 dF;
        dF.head<3>() = -b.m * dwwc;
        dF.tail<3>() = -(dw.cross(b.I * w) + w.cross(b.I * dw)) - b.m * b.c.cross(dwwc);
        return dF;
    }

    static Vec6 toEigen(const Deriv &d) {
        Vec6 v;
        v << d.getVCenter()[0], d.getVCenter()[1], d.getVCenter()[2], d.getVOrientation()[0],
            d.getVOrientation()[1], d.getVOrientation()[2];
        return v;
    }
    static void add(Deriv &d, const Vec6 &v, double scale = 1.0) {
        for (int k = 0; k < 3; ++k) {
            d.getVCenter()[k] += scale * v[k];
            d.getVOrientation()[k] += scale * v[3 + k];
        }
    }

    Vec3 gravity() const {
        const auto g = this->getContext()->getGravity();
        return Vec3(g[0], g[1], g[2]);
    }

    // ---- Mass interface
    // ------------------------------------------------------------------------------------------
    using Inherit::addMDx;
    void addMDx(const sofa::core::MechanicalParams *, DataVecDeriv &f, const DataVecDeriv &dx,
                SReal factor) override {
        if (!this->isComponentStateValid()) return;
        sofa::helper::WriteAccessor<DataVecDeriv> res = f;
        sofa::helper::ReadAccessor<DataVecDeriv> d = dx;
        const auto x = this->mstate->readPositions();
        const std::size_t n =
            std::min<std::size_t>(res.size(), std::min<std::size_t>(d.size(), x.size()));
        for (std::size_t i = 0; i < n; ++i)
            add(res[i], massMatrix(body(i, x[i])) * toEigen(d[i]), factor);
    }

    using Inherit::accFromF;
    void accFromF(const sofa::core::MechanicalParams *, DataVecDeriv &a,
                  const DataVecDeriv &f) override {
        if (!this->isComponentStateValid()) return;
        sofa::helper::WriteOnlyAccessor<DataVecDeriv> acc = a;
        sofa::helper::ReadAccessor<DataVecDeriv> force = f;
        const auto x = this->mstate->readPositions();
        const std::size_t n =
            std::min<std::size_t>(acc.size(), std::min<std::size_t>(force.size(), x.size()));
        for (std::size_t i = 0; i < n; ++i) {
            const Vec6 r = massMatrix(body(i, x[i])).ldlt().solve(toEigen(force[i]));
            acc[i] = Deriv();
            add(acc[i], r);
        }
    }

    /// Weight AND the velocity-dependent inertial wrench (gyroscopic, centrifugal): the exact
    /// Newton-Euler force.
    void addForce(const sofa::core::MechanicalParams *, DataVecDeriv &f, const DataVecCoord &x,
                  const DataVecDeriv &v) override {
        if (!this->isComponentStateValid()) return;
        sofa::helper::WriteAccessor<DataVecDeriv> force = f;
        sofa::helper::ReadAccessor<DataVecCoord> pos = x;
        sofa::helper::ReadAccessor<DataVecDeriv> vel = v;
        const Vec3 g = this->m_separateGravity.getValue() ? Vec3::Zero().eval() : gravity();
        const std::size_t n =
            std::min<std::size_t>(force.size(), std::min<std::size_t>(pos.size(), vel.size()));
        for (std::size_t i = 0; i < n; ++i) {
            const Vec6 vi = toEigen(vel[i]);
            add(force[i], wrench(body(i, pos[i]), vi.tail<3>(), g));
        }
    }

    void addDForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &df,
                   const DataVecDeriv &dx) override {
        if (!this->isComponentStateValid()) return;
        const double kFactor = mparams->kFactor(), bFactor = mparams->bFactor();
        const bool k = d_stiffness.getValue() && kFactor != 0.0,
                   b = d_gyroDamping.getValue() && bFactor != 0.0;
        if (!k && !b) return;
        sofa::helper::WriteAccessor<DataVecDeriv> res = df;
        sofa::helper::ReadAccessor<DataVecDeriv> d = dx;
        const auto x = this->mstate->readPositions();
        const auto v = this->mstate->readVelocities();
        const Vec3 g = this->m_separateGravity.getValue() ? Vec3::Zero().eval() : gravity();
        const std::size_t n = std::min<std::size_t>({res.size(), d.size(), x.size(), v.size()});
        for (std::size_t i = 0; i < n; ++i) {
            const Body bd = body(i, x[i]);
            const Vec6 di = toEigen(d[i]);
            const Vec3 w = toEigen(v[i]).tail<3>();
            if (k) add(res[i], wrenchRotationDerivative(bd, w, g, di.tail<3>()), kFactor);
            if (b) add(res[i], wrenchVelocityDerivative(bd, w, di.tail<3>()), bFactor);
        }
    }

    SReal getKineticEnergy(const sofa::core::MechanicalParams *,
                           const DataVecDeriv &v) const override {
        if (!this->isComponentStateValid()) return 0;
        sofa::helper::ReadAccessor<DataVecDeriv> vel = v;
        const auto x = this->mstate->readPositions();
        double e = 0;
        for (std::size_t i = 0; i < std::min<std::size_t>(vel.size(), x.size()); ++i) {
            const Vec6 vi = toEigen(vel[i]);
            e += 0.5 * vi.dot(massMatrix(body(i, x[i])) * vi);
        }
        return e;
    }

    SReal getPotentialEnergy(const sofa::core::MechanicalParams *,
                             const DataVecCoord &x) const override {
        if (!this->isComponentStateValid()) return 0;
        sofa::helper::ReadAccessor<DataVecCoord> pos = x;
        const Vec3 g = gravity();
        double e = 0;
        for (std::size_t i = 0; i < pos.size(); ++i) {
            const Body b = body(i, pos[i]);
            const Vec3 p(pos[i].getCenter()[0], pos[i].getCenter()[1], pos[i].getCenter()[2]);
            e -= b.m * g.dot(p + b.c);
        }
        return e;
    }

    sofa::type::Vec6 getMomentum(const sofa::core::MechanicalParams *, const DataVecCoord &x,
                                 const DataVecDeriv &v) const override {
        sofa::type::Vec6 out;
        if (!this->isComponentStateValid()) return out;
        sofa::helper::ReadAccessor<DataVecCoord> pos = x;
        sofa::helper::ReadAccessor<DataVecDeriv> vel = v;
        Vec6 total = Vec6::Zero();
        for (std::size_t i = 0; i < std::min<std::size_t>(pos.size(), vel.size()); ++i) {
            const Body b = body(i, pos[i]);
            const Vec6 vi = toEigen(vel[i]);
            const Vec3 vc =
                vi.head<3>() + vi.tail<3>().cross(b.c);  // velocity of the centre of mass
            const Vec3 pc =
                Vec3(pos[i].getCenter()[0], pos[i].getCenter()[1], pos[i].getCenter()[2]) + b.c;
            total.head<3>() += b.m * vc;
            total.tail<3>() += pc.cross(b.m * vc) + b.I * vi.tail<3>();
        }
        for (int k = 0; k < 6; ++k) out[k] = total[k];
        return out;
    }

    void addGravityToV(const sofa::core::MechanicalParams *mparams, DataVecDeriv &d_v) override {
        if (!this->isComponentStateValid() || this->m_separateGravity.getValue() == false) return;
        sofa::helper::WriteAccessor<DataVecDeriv> v = d_v;
        const Vec3 g = gravity() * mparams->dt();
        for (std::size_t i = 0; i < v.size(); ++i)
            for (int k = 0; k < 3; ++k)
                v[i].getVCenter()[k] +=
                    g[k];  // free fall: same acceleration for every point, no rotation
    }

    // ---- assembled matrices (direct solvers)
    // -----------------------------------------------------------------------
    using Inherit::addMToMatrix;
    void addMToMatrix(sofa::linearalgebra::BaseMatrix *matrix, SReal mFact,
                      unsigned int &offset) override {
        if (!this->isComponentStateValid()) return;
        const auto x = this->mstate->readPositions();
        for (std::size_t i = 0; i < x.size(); ++i) {
            const Mat6 M = massMatrix(body(i, x[i]));
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c)
                    matrix->add(offset + 6 * i + r, offset + 6 * i + c, mFact * M(r, c));
        }
    }

    void buildMassMatrix(sofa::core::behavior::MassMatrixAccumulator *matrices) override {
        if (!this->isComponentStateValid()) return;
        const auto x = this->mstate->readPositions();
        for (std::size_t i = 0; i < x.size(); ++i) {
            const Mat6 M = massMatrix(body(i, x[i]));
            for (int r = 0; r < 6; ++r)
                for (int c = 0; c < 6; ++c) matrices->add(6 * i + r, 6 * i + c, M(r, c));
        }
    }

    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix *matrix) override {
        if (!this->isComponentStateValid() || !d_stiffness.getValue()) return;
        auto dfdx =
            matrix->getForceDerivativeIn(this->mstate).withRespectToPositionsIn(this->mstate);
        const auto x = this->mstate->readPositions();
        const auto v = this->mstate->readVelocities();
        const Vec3 g = this->m_separateGravity.getValue() ? Vec3::Zero().eval() : gravity();
        for (std::size_t i = 0; i < std::min<std::size_t>(x.size(), v.size()); ++i) {
            const Body b = body(i, x[i]);
            const Vec3 w(v[i].getVOrientation()[0], v[i].getVOrientation()[1],
                         v[i].getVOrientation()[2]);
            for (int c = 0; c < 3; ++c) {  // columns: rotation increments
                const Vec6 col = wrenchRotationDerivative(b, w, g, Vec3::Unit(c));
                for (int r = 0; r < 6; ++r) dfdx(6 * i + r, 6 * i + 3 + c) += col[r];
            }
        }
    }

    void buildDampingMatrix(sofa::core::behavior::DampingMatrix *matrix) override {
        if (!this->isComponentStateValid() || !d_gyroDamping.getValue()) return;
        auto dfdv =
            matrix->getForceDerivativeIn(this->mstate).withRespectToVelocityIn(this->mstate);
        const auto x = this->mstate->readPositions();
        const auto v = this->mstate->readVelocities();
        for (std::size_t i = 0; i < std::min<std::size_t>(x.size(), v.size()); ++i) {
            const Body b = body(i, x[i]);
            const Vec3 w(v[i].getVOrientation()[0], v[i].getVOrientation()[1],
                         v[i].getVOrientation()[2]);
            for (int c = 0; c < 3; ++c) {
                const Vec6 col = wrenchVelocityDerivative(b, w, Vec3::Unit(c));
                for (int r = 0; r < 6; ++r) dfdv(6 * i + r, 6 * i + 3 + c) += col[r];
            }
        }
    }

    SReal getElementMass(sofa::Index index) const override {
        const auto &m = d_mass.getValue();
        return m.empty() ? 0.0 : m[m.size() == 1 ? 0 : std::min<std::size_t>(index, m.size() - 1)];
    }
    void getElementMass(sofa::Index index, sofa::linearalgebra::BaseMatrix *m) const override {
        m->resize(6, 6);
        const auto x = this->mstate->readPositions();
        if (index >= x.size()) return;
        const Mat6 M = massMatrix(body(index, x[index]));
        for (int r = 0; r < 6; ++r)
            for (int c = 0; c < 6; ++c) m->set(r, c, M(r, c));
    }
    bool isDiagonal() const override { return false; }

   private:
    std::size_t m_n = 0;
};

}  // namespace Cosserat::mapping
