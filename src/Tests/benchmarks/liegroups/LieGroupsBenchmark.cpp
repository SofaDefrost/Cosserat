#include <benchmark/benchmark.h>
#include <Cosserat/liegroups/RealSpace.h>
#include <Cosserat/liegroups/SO2.h>
#include <Cosserat/liegroups/SE2.h>
#include <Eigen/Core>

// RealSpace benchmarks
static void BM_RealSpace_Composition(benchmark::State& state) {
    // Setup
    const int dim = state.range(0);
    Cosserat::RealSpace<double, Eigen::Dynamic> a(Eigen::VectorXd::Random(dim));
    Cosserat::RealSpace<double, Eigen::Dynamic> b(Eigen::VectorXd::Random(dim));
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.compose(b));
    }
    
    state.SetComplexityN(dim);
}
BENCHMARK(BM_RealSpace_Composition)->RangeMultiplier(2)->Range(1, 1024)->Complexity();

static void BM_RealSpace_Inverse(benchmark::State& state) {
    // Setup
    const int dim = state.range(0);
    Cosserat::RealSpace<double, Eigen::Dynamic> a(Eigen::VectorXd::Random(dim));
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.inverse());
    }
    
    state.SetComplexityN(dim);
}
BENCHMARK(BM_RealSpace_Inverse)->RangeMultiplier(2)->Range(1, 1024)->Complexity();

// SO2 benchmarks
static void BM_SO2_Composition(benchmark::State& state) {
    // Setup
    Cosserat::SO2<double> a(M_PI / 4.0);
    Cosserat::SO2<double> b(M_PI / 3.0);
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.compose(b));
    }
}
BENCHMARK(BM_SO2_Composition);

static void BM_SO2_Inverse(benchmark::State& state) {
    // Setup
    Cosserat::SO2<double> a(M_PI / 4.0);
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.inverse());
    }
}
BENCHMARK(BM_SO2_Inverse);

static void BM_SO2_Log(benchmark::State& state) {
    // Setup
    Cosserat::SO2<double> a(M_PI / 4.0);
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.log());
    }
}
BENCHMARK(BM_SO2_Log);

static void BM_SO2_Exp(benchmark::State& state) {
    // Setup
    double theta = M_PI / 4.0;
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(Cosserat::SO2<double>::exp(theta));
    }
}
BENCHMARK(BM_SO2_Exp);

// SE2 benchmarks
static void BM_SE2_Composition(benchmark::State& state) {
    // Setup
    Cosserat::SE2<double> a(M_PI / 4.0, 1.0, 2.0);
    Cosserat::SE2<double> b(M_PI / 3.0, -1.0, 0.5);
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.compose(b));
    }
}
BENCHMARK(BM_SE2_Composition);

static void BM_SE2_Inverse(benchmark::State& state) {
    // Setup
    Cosserat::SE2<double> a(M_PI / 4.0, 1.0, 2.0);
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.inverse());
    }
}
BENCHMARK(BM_SE2_Inverse);

static void BM_SE2_Log(benchmark::State& state) {
    // Setup
    Cosserat::SE2<double> a(M_PI / 4.0, 1.0, 2.0);
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(a.log());
    }
}
BENCHMARK(BM_SE2_Log);

static void BM_SE2_Exp(benchmark::State& state) {
    // Setup
    Eigen::Vector3d tangent;
    tangent << M_PI / 4.0, 1.0, 2.0;
    
    // Benchmark loop
    for (auto _ : state) {
        benchmark::DoNotOptimize(Cosserat::SE2<double>::exp(tangent));
    }
}
BENCHMARK(BM_SE2_Exp);

BENCHMARK_MAIN();

