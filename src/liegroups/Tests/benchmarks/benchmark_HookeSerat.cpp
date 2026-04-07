#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <benchmark/benchmark.h>
#include <liegroups/SE3.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

using namespace sofa::component::cosserat::liegroups;
using namespace Cosserat::mapping;

class ConcreteHookeSeratMapping
	: public HookeSeratBaseMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
								   sofa::defaulttype::Rigid3Types> {
public:
	void doBaseCosseratInit() override {}
	using HookeSeratBaseMapping::computeTangExpImplementation;
};

static void BM_JacobianComputation(benchmark::State &state) {
	ConcreteHookeSeratMapping::TangentVector strain;
	strain << 0.1, 0.2, 0.3, 1.0, 0.1, 0.0;
	double curv_abs = 1.0;

	ConcreteHookeSeratMapping::AdjointMatrix adjoint = ConcreteHookeSeratMapping::AdjointMatrix::Identity();
	ConcreteHookeSeratMapping::AdjointMatrix result;

	for (auto _: state) {
		ConcreteHookeSeratMapping::computeTangExpImplementation(curv_abs, strain, adjoint, result);
	}
}
BENCHMARK(BM_JacobianComputation);

static void BM_TrajectoryGeneration(benchmark::State &state) {
	ConcreteHookeSeratMapping mapping;
	SectionInfo section;
	section.setLength(1.0);
	ConcreteHookeSeratMapping::TangentVector strain;
	strain << 0.1, 0.0, 0.0, 1.0, 0.0, 0.0;
	section.setStrainsVec(strain);
	mapping.addSection(section);

	int num_points = state.range(0);

	for (auto _: state) {
		auto traj = mapping.generateSectionTrajectory(num_points);
		benchmark::DoNotOptimize(traj);
	}
}
BENCHMARK(BM_TrajectoryGeneration)->Range(10, 1000);

BENCHMARK_MAIN();
