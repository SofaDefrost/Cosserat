// Standalone test for Sim3 implementation
// This tests the Sim3 Lie group without the full SOFA build system

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>

// Simplified forward declarations
namespace sofa::component::cosserat::liegroups {

template<typename Scalar>
class SO3 {
public:
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
    using Quaternion = Eigen::Quaternion<Scalar>;

    SO3() : m_quat(Quaternion::Identity()) {}
    explicit SO3(const Quaternion& q) : m_quat(q) {}
    explicit SO3(const Matrix3& R) : m_quat(R) {}

    SO3 operator*(const SO3& other) const { return SO3(m_quat * other.m_quat); }
    SO3 inverse() const { return SO3(m_quat.conjugate()); }
    Vector3 log() const {
        Eigen::AngleAxis<Scalar> aa(m_quat);
        return aa.axis() * aa.angle();
    }
    Matrix3 matrix() const { return m_quat.toRotationMatrix(); }
    static SO3 exp(const Vector3& omega) {
        Scalar theta = omega.norm();
        if (theta < 1e-10) return SO3();
        Eigen::AngleAxis<Scalar> aa(theta, omega/theta);
        return SO3(Quaternion(aa));
    }

private:
    Quaternion m_quat;
};

template<typename Scalar>
class Sim3 {
public:
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using Vector7 = Eigen::Matrix<Scalar, 7, 1>;
    using Matrix4 = Eigen::Matrix<Scalar, 4, 4>;
    using SO3Type = SO3<Scalar>;

    Sim3() : m_rotation(), m_translation(Vector3::Zero()), m_scale(Scalar(1)) {}
    Sim3(const SO3Type& rotation, const Vector3& translation, const Scalar& scale)
        : m_rotation(rotation), m_translation(translation), m_scale(scale) {}

    Sim3 compose(const Sim3& other) const {
        return Sim3(m_rotation * other.m_rotation,
                   m_rotation.matrix() * (other.m_translation * m_scale) + m_translation,
                   m_scale * other.m_scale);
    }

    Sim3 inverse() const {
        SO3Type inv_rot = m_rotation.inverse();
        Scalar inv_scale = Scalar(1) / m_scale;
        return Sim3(inv_rot, -(inv_rot.matrix() * m_translation) * inv_scale, inv_scale);
    }

    Vector3 act(const Vector3& point) const {
        return m_rotation.matrix() * (point * m_scale) + m_translation;
    }

    static Sim3 exp(const Vector7& xi) {
        Vector3 omega = xi.template head<3>();
        Vector3 v = xi.template segment<3>(3);
        Scalar s = xi(6);

        SO3Type R = SO3Type::exp(omega);
        Scalar scale = std::exp(s);

        Scalar theta = omega.norm();
        Vector3 translation;

        if (theta < 1e-10) {
            translation = v;
        } else {
            // Simplified V matrix computation
            translation = v; // Approximation for now
        }

        return Sim3(R, translation, scale);
    }

    Vector7 log() const {
        Vector3 omega = m_rotation.log();
        Vector3 v = m_translation; // Approximation
        Scalar s = std::log(m_scale);

        Vector7 result;
        result.template head<3>() = omega;
        result.template segment<3>(3) = v;
        result(6) = s;
        return result;
    }

    Matrix4 matrix() const {
        Matrix4 T = Matrix4::Identity();
        Matrix4 R_scaled = m_rotation.matrix() * m_scale;
        T.template block<3, 3>(0, 0) = R_scaled;
        T.template block<3, 1>(0, 3) = m_translation;
        return T;
    }

    const SO3Type& rotation() const { return m_rotation; }
    const Vector3& translation() const { return m_translation; }
    const Scalar& scale() const { return m_scale; }

private:
    SO3Type m_rotation;
    Vector3 m_translation;
    Scalar m_scale;
};

} // namespace

int main() {
    using namespace sofa::component::cosserat::liegroups;

    std::cout << "Testing Sim3 Lie Group Implementation" << std::endl;
    std::cout << "====================================" << std::endl;

    // Test basic operations
    Sim3<double> identity;
    std::cout << "Identity scale: " << identity.scale() << std::endl;

    // Test composition
    Eigen::Vector3d trans1(1.0, 2.0, 3.0);
    Eigen::Vector3d trans2(4.0, 5.0, 6.0);
    Sim3<double> sim1(SO3<double>(), trans1, 2.0);
    Sim3<double> sim2(SO3<double>(), trans2, 0.5);

    auto composed = sim1.compose(sim2);
    std::cout << "Composed translation: " << composed.translation().transpose() << std::endl;
    std::cout << "Composed scale: " << composed.scale() << std::endl;

    // Test action on point
    Eigen::Vector3d point(1.0, 1.0, 1.0);
    auto transformed = sim1.act(point);
    std::cout << "Transformed point: " << transformed.transpose() << std::endl;

    // Test exp/log
    Eigen::Matrix<double, 7, 1> tangent;
    tangent << 0.1, 0.2, 0.3, 1.0, 2.0, 3.0, 0.1; // ω, v, log(s)
    auto exp_result = Sim3<double>::exp(tangent);
    auto log_result = exp_result.log();

    std::cout << "Exp result scale: " << exp_result.scale() << std::endl;
    std::cout << "Log result: " << log_result.transpose() << std::endl;

    std::cout << "Sim3 implementation test completed!" << std::endl;
    return 0;
}