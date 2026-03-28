#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.inl>
#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

using namespace Cosserat::mapping;
using namespace sofa::component::cosserat::liegroups;

// Concrete subclass for testing
class ConcreteHookeSeratMapping
	: public HookeSeratBaseMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
								   sofa::defaulttype::Rigid3Types> {
public:
	using In1 = sofa::defaulttype::Vec3Types;
	using In2 = sofa::defaulttype::Rigid3Types;
	using Out = sofa::defaulttype::Rigid3Types;

	void doBaseCosseratInit() override {}

	// Implement pure virtual methods from Multi2Mapping
	void applyDJT(const sofa::core::MechanicalParams * /*mparams*/, sofa::core::MultiVecDerivId /*inForce*/,
				  sofa::core::ConstMultiVecDerivId /*outForce*/) override {}

	void apply(const sofa::core::MechanicalParams * /*mparams*/,
			   const sofa::type::vector<sofa::DataVecCoord_t<Out> *> & /*dataVecOutPos*/,
			   const sofa::type::vector<const sofa::DataVecCoord_t<In1> *> & /*dataVecIn1Pos*/,
			   const sofa::type::vector<const sofa::DataVecCoord_t<In2> *> & /*dataVecIn2Pos*/) override {}

	void applyJT(const sofa::core::MechanicalParams * /*mparams*/,
				 const sofa::type::vector<sofa::DataVecDeriv_t<In1> *> & /*dataVecOut1Force*/,
				 const sofa::type::vector<sofa::DataVecDeriv_t<In2> *> & /*dataVecOut2RootForce*/,
				 const sofa::type::vector<const sofa::DataVecDeriv_t<Out> *> & /*dataVecInForce*/) override {}
};

class HookeSeratPhase1Test : public ::testing::Test {
protected:
	using SE3Type = SE3<double>;
	using TangentVector = typename SE3Type::TangentVector;
	using Vector3 = typename SE3Type::Vector3;

	ConcreteHookeSeratMapping mapping;

	void SetUp() override {
		// Clear any existing data
		mapping.clearSections();
	}
};

TEST_F(HookeSeratPhase1Test, ValidateJacobianAccuracy) {
	// Add a section with some strain
	TangentVector strain;
	strain << 0.1, 0.05, 0.02, 0.01, 0.01, 0.01;
	SectionInfo section(1.0, strain, 0);
	mapping.addSection(section);

	// Validate
	EXPECT_TRUE(mapping.validateJacobianAccuracy(1e-4));
}

TEST_F(HookeSeratPhase1Test, GenerateSectionTrajectory) {
	// Add a section
	TangentVector strain;
	strain << 0.1, 0.0, 0.0, 1.0, 0.0, 0.0; // Constant curvature and elongation
	SectionInfo section(1.0, strain, 0);
	mapping.addSection(section);

	int num_points = 5;
	auto trajectory = mapping.generateSectionTrajectory(num_points);

	// Check size: num_points * sections + 1 (start to end)
	// Here 1 section, so 5 + 1 = 6 points
	EXPECT_EQ(trajectory.size(), num_points + 1);

	// Check first point (t=0)
	// The transformation should be Identity (assuming section starts at Identity)
	EXPECT_TRUE(trajectory[0].getTransformation().isApprox(SE3Type::computeIdentity()));

	// Check last point (t=1)
	SE3Type expected_end = section.getLocalTransformation(1.0);
	EXPECT_TRUE(trajectory.back().getTransformation().isApprox(expected_end));

	// Check intermediate point
	// t = 0.5 (index 2 if num_points=4? No, num_points=5.
	// Indices: 0(0.0), 1(0.2), 2(0.4), 3(0.6), 4(0.8), 5(1.0) ?
	// Loop goes 0 to num_points-1.
	// i=0 -> t=0.0
	// i=1 -> t=0.2
	// ...
	// i=4 -> t=0.8
	// Then we add last point t=1.0.
	// So yes.

	// Check index 2 (t=0.4)
	double t = 0.4;
	SE3Type expected_mid = section.getLocalTransformation(t);
	EXPECT_TRUE(trajectory[2].getTransformation().isApprox(expected_mid));
}
