/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                  *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* along with this program. If not, see <http://www.gnu.org/licenses/\>.        *
******************************************************************************/

#include <sofa/testing/BaseTest.h>
#include <benchmark/benchmark.h>
#include <Cosserat/liegroups/Bundle.h>
#include <Cosserat/liegroups/RealSpace.h>
#include <Cosserat/liegroups/SO3.h>
#include <Cosserat/liegroups/SE3.h>
#include <Cosserat/liegroups/SE23.h>
#include <Cosserat/liegroups/SGal3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <random>
#include <vector>

namespace sofa::component::cosserat::liegroups::benchmark {

using Vector3 = Eigen::Vector3d;
using Matrix3 = Eigen::Matrix3d;
using Quaternion = Eigen::Quaterniond;

// Helper for random generation
class RandomGenerator {
private:
    std::mt19937 gen{std::random_device{}()};
    std::uniform_real_distribution<double> angle_dist{0, 2*M_PI};
    std::normal_distribution<double> vec_dist{0, 1.0};

public:
    Vector3 randomVector() {
        return Vector3(vec_dist(gen), vec_dist(gen), vec_dist(gen));
    }

    Vector3 randomUnitVector() {
        Vector3 v = randomVector();
        return v.normalized();
    }

    double randomAngle() {
        return angle_dist(gen);
    }
};

/**
 * Benchmark SO3 operations
 */
static void BM_SO3_Operations(benchmark::State& state) {
    RandomGenerator rng;
    SO3<double> rot(rng.randomAngle(), rng.randomUnitVector());
    Vector3 point = rng.randomVector();

    for (auto _ : state) {
        // Test common operations
        auto result1 = rot.act(point);
        auto result2 = rot.inverse();
        auto result3 = rot.log();
        auto result4 = rot.adjoint();
        benchmark::DoNotOptimize(result1);
        benchmark::DoNotOptimize(result2);
        benchmark::DoNotOptimize(result3);
        benchmark::DoNotOptimize(result4);
    }
}
BENCHMARK(BM_SO3_Operations);

/**
 * Benchmark SE3 operations
 */
static void BM_SE3_Operations(benchmark::State& state) {
    RandomGenerator rng;
    SE3<double> transform(
        SO3<double>(rng.randomAngle(), rng.randomUnitVector()),
        rng.randomVector()
    );
    Vector3 point = rng.randomVector();

    for (auto _ : state) {
        // Test common operations
        auto result1 = transform.act(point);
        auto result2 = transform.inverse();
        auto result3 = transform.log();
        auto result4 = transform.adjoint();
        benchmark::DoNotOptimize(result1);
        benchmark::DoNotOptimize(result2);
        benchmark::DoNotOptimize(result3);
        benchmark::DoNotOptimize(result4);
    }
}
BENCHMARK(BM_SE3_Operations);

/**
 * Benchmark SE_2(3) operations
 */
static void BM_SE23_Operations(benchmark::State& state) {
    RandomGenerator rng;
    SE23<double> extended_pose(
        SE3<double>(
            SO3<double>(rng.randomAngle(), rng.randomUnitVector()),
            rng.randomVector()
        ),
        rng.randomVector()
    );
    Vector3 point = rng.randomVector();

    for (auto _ : state) {
        // Test common operations
        auto result1 = extended_pose.act(point);
        auto result2 = extended_pose.inverse();
        auto result3 = extended_pose.log();
        auto result4 = extended_pose.adjoint();
        benchmark::DoNotOptimize(result1);
        benchmark::DoNotOptimize(result2);
        benchmark::DoNotOptimize(result3);
        benchmark::DoNotOptimize(result4);
    }
}
BENCHMARK(BM_SE23_Operations);

/**
 * Benchmark Bundle operations
 */
static void BM_Bundle_Operations(benchmark::State& state) {
    RandomGenerator rng;
    using PoseVel = Bundle<SE3<double>, RealSpace<double, 3>>;
    
    PoseVel bundle(
        SE3<double>(
            SO3<double>(rng.randomAngle(), rng.randomUnitVector()),
            rng.randomVector()
        ),
        RealSpace<double, 3>(rng.randomVector())
    );

    for (auto _ : state) {
        // Test common operations
        auto result1 = bundle.inverse();
        auto result2 = bundle.log();
        auto result3 = bundle.adjoint();
        benchmark::DoNotOptimize(result1);
        benchmark::DoNotOptimize(result2);
        benchmark::DoNotOptimize(result3);
    }
}
BENCHMARK(BM_Bundle_Operations);

/**
 * Benchmark Cosserat rod operations
 */
static void BM_CosseratRod_Operations(benchmark::State& state) {
    RandomGenerator rng;
    const int num_segments = state.range(0);
    using RodSegment = Bundle<SE3<double>, RealSpace<double, 3>>;
    std::vector<RodSegment> segments;

    // Initialize rod segments
    for (int i = 0; i < num_segments; ++i) {
        segments.push_back(RodSegment(
            SE3<double>(
                SO3<double>(rng.randomAngle(), rng.randomUnitVector()),
                rng.randomVector()
            ),
            RealSpace<double, 3>(rng.randomVector())
        ));
    }

    for (auto _ : state) {
        // Simulate rod deformation
        for (int i = 1; i < num_segments; ++i) {
            auto rel_transform = segments[i-1].inverse() * segments[i];
            auto strain = rel_transform.log();
            benchmark::DoNotOptimize(strain);
        }
    }
}
BENCHMARK(BM_CosseratRod_Operations)
    ->RangeMultiplier(2)
    ->Range(8, 128);

/**
 * Benchmark exponential map implementations
 */
static void BM_ExpMap_Operations(benchmark::State& state) {
    RandomGenerator rng;
    Vector3 omega = rng.randomVector();

    for (auto _ : state) {
        // SO3 exponential
        auto rot = SO3<double>().exp(omega);
        
        // SE3 exponential
        Vector3 v = rng.randomVector();
        auto transform = SE3<double>().exp(
            (Eigen::Matrix<double, 6, 1>() << v, omega).finished()
        );
        
        benchmark::DoNotOptimize(rot);
        benchmark::DoNotOptimize(transform);
    }
}
BENCHMARK(BM_ExpMap_Operations);

/**
 * Benchmark interpolation operations
 */
static void BM_Interpolation_Operations(benchmark::State& state) {
    RandomGenerator rng;
    
    // Create random transformations
    SE3<double> T1(
        SO3<double>(rng.randomAngle(), rng.randomUnitVector()),
        rng.randomVector()
    );
    SE3<double> T2(
        SO3<double>(rng.randomAngle(), rng.randomUnitVector()),
        rng.randomVector()
    );

    const int num_steps = state.range(0);
    std::vector<double> times(num_steps);
    for (int i = 0; i < num_steps; ++i) {
        times[i] = static_cast<double>(i) / (num_steps - 1);
    }

    for (auto _ : state) {
        // Interpolate between transformations
        for (double t : times) {
            auto result = interpolate(T1, T2, t);
            benchmark::DoNotOptimize(result);
        }
    }
}
BENCHMARK(BM_Interpolation_Operations)
    ->RangeMultiplier(2)
    ->Range(8, 128);

} // namespace sofa::component::cosserat::liegroups::benchmark

BENCHMARK_MAIN();
