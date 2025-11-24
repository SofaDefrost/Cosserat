/**
 * Test program to compare the new Lie group implementation vs legacy trigonometric implementation
 * of the tangent exponential map in HookeSeratBaseMapping.
 *
 * This test verifies that SE3Type::rightJacobian() produces equivalent results to the
 * manual trigonometric series expansion for various strain configurations.
 */

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <liegroups/SE3.h>
#include <liegroups/Types.h>

namespace sofa::component::cosserat::liegroups {
    using SE3Type = SE3<double>;
    using SO3Type = SO3<double>;
    using Vector3 = typename SE3Type::Vector3;
    using TangentVector = typename SE3Type::TangentVector;
    using AdjointMatrix = typename SE3Type::AdjointMatrix;
    using Matrix3 = typename SE3Type::Matrix3;
    using Matrix4 = typename SE3Type::Matrix4;
}

/**
 * Legacy implementation using manual trigonometric series
 */
void computeTangExpImplementationLegacy(const double& curv_abs,
    const sofa::component::cosserat::liegroups::TangentVector & strain,
    const sofa::component::cosserat::liegroups::AdjointMatrix &adjoint_matrix,
    sofa::component::cosserat::liegroups::AdjointMatrix & tang_adjoint_matrix)
{
    using namespace sofa::component::cosserat::liegroups;
    using SReal = double;

    SReal theta = Vector3(strain(0), strain(1), strain(2)).norm();

    tang_adjoint_matrix = AdjointMatrix::Zero();
    AdjointMatrix Id6 = AdjointMatrix::Identity();

    if (theta <= std::numeric_limits<double>::epsilon()) {
        double scalar0 = std::pow(curv_abs, 2) / 2.0;
        tang_adjoint_matrix = curv_abs * Id6 + scalar0 * adjoint_matrix;
    } else {
        double scalar1 = (4.0 - 4.0 * cos(curv_abs * theta) -
                          curv_abs * theta * sin(curv_abs * theta)) /
                         (2.0 * theta * theta);
        double scalar2 = (4.0 * curv_abs * theta +
                          curv_abs * theta * cos(curv_abs * theta) -
                          5.0 * sin(curv_abs * theta)) /
                         (2.0 * theta * theta * theta);
        double scalar3 = (2.0 - 2.0 * cos(curv_abs * theta) -
                          curv_abs * theta * sin(curv_abs * theta)) /
                         (2.0 * theta * theta * theta * theta);
        double scalar4 = (2.0 * curv_abs * theta +
                          curv_abs * theta * cos(curv_abs * theta) -
                          3.0 * sin(curv_abs * theta)) /
                         (2.0 * theta * theta * theta * theta * theta);

        tang_adjoint_matrix = curv_abs * Id6 + scalar1 * adjoint_matrix + scalar2 * adjoint_matrix * adjoint_matrix +
              scalar3 * adjoint_matrix * adjoint_matrix * adjoint_matrix +
              scalar4 * adjoint_matrix * adjoint_matrix * adjoint_matrix * adjoint_matrix;
    }
}

/**
 * New implementation using Lie group right Jacobian
 */
void computeTangExpImplementationNew(const double& curv_abs,
    const sofa::component::cosserat::liegroups::TangentVector & strain,
    const sofa::component::cosserat::liegroups::AdjointMatrix & /*adjoint_matrix*/,
    sofa::component::cosserat::liegroups::AdjointMatrix & tang_adjoint_matrix)
{
    using namespace sofa::component::cosserat::liegroups;

    TangentVector scaled_strain = strain * curv_abs;
    tang_adjoint_matrix = SE3Type::rightJacobian(scaled_strain);
}

/**
 * Test function to compare implementations
 */
bool testEquivalence(const double& curv_abs,
    const sofa::component::cosserat::liegroups::TangentVector & strain,
    const sofa::component::cosserat::liegroups::AdjointMatrix &adjoint_matrix,
    double tolerance = 1e-6)
{
    using namespace sofa::component::cosserat::liegroups;

    AdjointMatrix new_result, legacy_result;

    // Compute using new Lie group implementation
    computeTangExpImplementationNew(curv_abs, strain, adjoint_matrix, new_result);

    // Compute using legacy trigonometric implementation
    computeTangExpImplementationLegacy(curv_abs, strain, adjoint_matrix, legacy_result);

    // Compare results
    AdjointMatrix diff = new_result - legacy_result;
    double max_diff = diff.cwiseAbs().maxCoeff();

    std::cout << "Strain: [" << strain.transpose() << "]" << std::endl;
    std::cout << "Curvilinear abs: " << curv_abs << std::endl;
    std::cout << "Max difference: " << max_diff << std::endl;

    if (max_diff > tolerance) {
        std::cout << "❌ FAILED - Difference exceeds tolerance " << tolerance << std::endl;
        std::cout << "New result:\n" << new_result << std::endl;
        std::cout << "Legacy result:\n" << legacy_result << std::endl;
        std::cout << "Difference:\n" << diff << std::endl;
        return false;
    }

    std::cout << "✅ PASSED" << std::endl << std::endl;
    return true;
}

int main() {
    using namespace sofa::component::cosserat::liegroups;

    std::cout << "Testing equivalence between new Lie group and legacy trigonometric implementations" << std::endl;
    std::cout << "=================================================================================" << std::endl;

    // Test cases with different strain configurations
    std::vector<std::pair<double, TangentVector>> test_cases = {
        // Small strain case
        {0.1, TangentVector{0.01, 0.01, 0.01, 0.0, 0.0, 0.0}},

        // Medium strain case
        {0.5, TangentVector{0.1, 0.2, 0.15, 0.05, 0.03, 0.08}},

        // Large strain case
        {1.0, TangentVector{0.5, 0.3, 0.4, 0.2, 0.1, 0.15}},

        // Pure rotation case
        {0.8, TangentVector{0.3, 0.2, 0.1, 0.0, 0.0, 0.0}},

        // Pure translation case
        {0.6, TangentVector{0.0, 0.0, 0.0, 0.4, 0.3, 0.2}},

        // Zero strain case
        {0.2, TangentVector{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}},
    };

    // Sample adjoint matrix (identity for simplicity, but could be any valid adjoint)
    AdjointMatrix adjoint_matrix = AdjointMatrix::Identity();

    bool all_passed = true;
    for (const auto& [curv_abs, strain] : test_cases) {
        bool passed = testEquivalence(curv_abs, strain, adjoint_matrix);
        all_passed &= passed;
    }

    std::cout << "=================================================================================" << std::endl;
    if (all_passed) {
        std::cout << "🎉 ALL TESTS PASSED - Implementations are equivalent!" << std::endl;
        return 0;
    } else {
        std::cout << "❌ SOME TESTS FAILED - Review implementation differences" << std::endl;
        return 1;
    }
}