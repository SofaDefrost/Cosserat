#include <Cosserat/mapping/CosseratGeometryMapping.h>
#include <benchmark/benchmark.h>
#include <liegroups/SE3.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

using namespace sofa::component::cosserat::liegroups;
using namespace Cosserat::mapping;

class ConcreteStrain2RigidCosseratMapping
	: public CosseratGeometryMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
								   sofa::defaulttype::Rigid3Types> {
public:
	void doBaseCosseratInit() override {}
	using CosseratGeometryMapping::computeTangExpImplementation;
};

static void BM_JacobianComputation(benchmark::State &state) {
	ConcreteStrain2RigidCosseratMapping::TangentVector strain;
	strain << 0.1, 0.2, 0.3, 1.0, 0.1, 0.0;
	double curv_abs = 1.0;

	ConcreteStrain2RigidCosseratMapping::AdjointMatrix result;

	for (auto _: state) {
		ConcreteStrain2RigidCosseratMapping::computeTangExpImplementation(curv_abs, strain, result);
	}
}
BENCHMARK(BM_JacobianComputation);

static void BM_TrajectoryGeneration(benchmark::State &state) {
	ConcreteStrain2RigidCosseratMapping mapping;
	SectionInfo section;
	section.setLength(1.0);
	ConcreteStrain2RigidCosseratMapping::TangentVector strain;
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
