#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <gtest/gtest.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

using namespace sofa::component::cosserat::liegroups;
using namespace Cosserat::mapping;

class ConcreteHookeSeratMapping
	: public HookeSeratBaseMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
								   sofa::defaulttype::Rigid3Types> {
public:
	using In1 = sofa::defaulttype::Vec3Types;
	using In2 = sofa::defaulttype::Rigid3Types;
	using Out = sofa::defaulttype::Rigid3Types;

	void doBaseCosseratInit() override {}

	// Implement pure virtual methods from Multi2Mapping
	void apply(const sofa::core::MechanicalParams * /* mparams */,
			   const sofa::type::vector<sofa::DataVecCoord_t<Out> *> & /* dataVecOutPos */,
			   const sofa::type::vector<const sofa::DataVecCoord_t<In1> *> & /* dataVecIn1Pos */,
			   const sofa::type::vector<const sofa::DataVecCoord_t<In2> *> & /* dataVecIn2Pos */) override {}

	void applyJ(const sofa::core::MechanicalParams * /* mparams */,
				const sofa::type::vector<sofa::DataVecDeriv_t<Out> *> & /* dataVecOutVel */,
				const sofa::type::vector<const sofa::DataVecDeriv_t<In1> *> & /* dataVecIn1Vel */,
				const sofa::type::vector<const sofa::DataVecDeriv_t<In2> *> & /* dataVecIn2Vel */) override {}

	void applyJT(const sofa::core::MechanicalParams * /* mparams */,
				 const sofa::type::vector<sofa::DataVecDeriv_t<In1> *> & /* dataVecOut1Force */,
				 const sofa::type::vector<sofa::DataVecDeriv_t<In2> *> & /* dataVecOut2Force */,
				 const sofa::type::vector<const sofa::DataVecDeriv_t<Out> *> & /* dataVecInForce */) override {}

	void applyJT(const sofa::core::ConstraintParams * /* cparams */,
				 const sofa::type::vector<sofa::DataMatrixDeriv_t<In1> *> & /* dataMatOut1Const */,
				 const sofa::type::vector<sofa::DataMatrixDeriv_t<In2> *> & /* dataMatOut2Const */,
				 const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out> *> & /* dataMatInConst */) override {}

	void applyDJT(const sofa::core::MechanicalParams * /* mparams */, sofa::core::MultiVecDerivId /* inForce */,
				  sofa::core::ConstMultiVecDerivId /* outForce */) override {}

	// Expose protected methods
	using HookeSeratBaseMapping::checkContinuity;
};

TEST(MultiSectionBeamTest, TopologyValidation) {
	ConcreteHookeSeratMapping mapping;

	// Invalid topology (size mismatch)
	BeamTopology invalid_topology;
	invalid_topology.parent_indices = {0};
	invalid_topology.relative_transforms.clear(); // Empty

	EXPECT_FALSE(invalid_topology.isValid());

	// Valid topology
	BeamTopology valid_topology;
	valid_topology.parent_indices = {-1, 0};
	valid_topology.relative_transforms.resize(2);
	valid_topology.connection_stiffnesses.resize(2);

	EXPECT_TRUE(valid_topology.isValid());

	mapping.setBeamTopology(valid_topology);
	EXPECT_EQ(mapping.getBeamTopology().getNumSections(), 2);
}

TEST(MultiSectionBeamTest, ContinuityCheck) {
	ConcreteHookeSeratMapping mapping;

	// Create two sections that are continuous
	// Section 1: Length 1, Strain 0 (Identity transform) -> End at Identity * Length? No, Exp(0)*L
	// If strain is 0, Exp(0) is Identity.
	// Wait, getLocalTransformation(t) = Exp(strain * s).
	// If strain is 0, transform is Identity for all t.

	SectionInfo s1;
	s1.setLength(1.0);
	s1.setStrain(TangentVector::Zero());

	SectionInfo s2;
	s2.setLength(1.0);
	s2.setStrain(TangentVector::Zero());

	mapping.addSection(s1);
	mapping.addSection(s2);

	// Check continuity
	// End of s1: Exp(0 * 1.0) = Identity
	// Start of s2: Exp(0 * 0.0) = Identity
	// Should be continuous
	EXPECT_TRUE(mapping.checkContinuity());

	// Create discontinuity
	SectionInfo s3;
	s3.setLength(1.0);
	TangentVector strain;
	strain << 1, 0, 0, 0, 0, 0;
	s3.setStrain(strain); // Non-zero strain

	mapping.clearSections();
	mapping.addSection(s3);
	mapping.addSection(s2); // s2 starts at Identity

	// End of s3: Exp(strain * 1.0) != Identity
	// Start of s2: Identity
	// Should be discontinuous
	EXPECT_FALSE(mapping.checkContinuity());
}
