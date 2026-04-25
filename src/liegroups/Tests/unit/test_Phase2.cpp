#include <Cosserat/mapping/BeamStateEstimator.h>
#include <Cosserat/mapping/CosseratGeometryMapping.h>
#include <gtest/gtest.h>
#include <liegroups/GaussianOnManifold.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

using namespace sofa::component::cosserat::liegroups;
using namespace Cosserat::mapping;

// Test GaussianOnManifold
TEST(GaussianOnManifoldTest, Initialization) {
	using SE3Type = SE3<double>;

	GaussianOnManifold<SE3Type> gaussian;

	EXPECT_TRUE(gaussian.getMean().isApprox(SE3Type::computeIdentity()));
	EXPECT_TRUE(gaussian.getCovariance().isApprox(Eigen::Matrix<double, 6, 6>::Identity()));
}

TEST(GaussianOnManifoldTest, Transformation) {
	using SE3Type = SE3<double>;
	using TangentVector = typename SE3Type::TangentVector;

	TangentVector v;
	v << 1, 0, 0, 0, 0, 0;
	SE3Type transform = SE3Type::computeExp(v);

	GaussianOnManifold<SE3Type> gaussian;
	auto transformed = gaussian.transform(transform);

	EXPECT_TRUE(transformed.getMean().isApprox(transform));
	// Covariance should be Ad(T) * I * Ad(T)^T
	// For pure translation, Ad(T) is identity for rotation part, but has off-diagonal for translation
	// Actually Ad(T) = [R, [t]xR; 0, R]
	// Here R=I, t=[1,0,0]
}

// Test BeamStateEstimator
TEST(BeamStateEstimatorTest, Prediction) {
	BeamStateEstimator estimator;
	BeamStateEstimator::SE3Type initial_pose = BeamStateEstimator::SE3Type::computeIdentity();
	BeamStateEstimator::CovarianceMatrix initial_cov = BeamStateEstimator::CovarianceMatrix::Identity() * 0.1;

	estimator.initialize(initial_pose, initial_cov);

	BeamStateEstimator::TangentVector strain;
	strain << 0.1, 0.0, 0.0, 1.0, 0.0, 0.0; // Constant curvature and elongation
	double dt = 1.0;
	BeamStateEstimator::CovarianceMatrix process_noise = BeamStateEstimator::CovarianceMatrix::Identity() * 0.01;

	estimator.predict(strain, dt, process_noise);

	// Mean should have moved
	EXPECT_FALSE(estimator.getEstimate().getMean().isApprox(initial_pose));

	// Uncertainty should have increased
	EXPECT_GT(estimator.getEstimationConfidence(), initial_cov.trace());
}

// Test BeamTopology
class ConcreteStrain2RigidCosseratMapping
	: public CosseratGeometryMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
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
};

TEST(BeamTopologyTest, Structure) {
	ConcreteStrain2RigidCosseratMapping mapping;

	// Create a Y-shape topology: 0 -> 1, 0 -> 2
	BeamTopology topology;
	topology.parent_indices = {-1, 0, 0}; // 0 is root, 1 and 2 are children of 0
	topology.relative_transforms.resize(3);
	topology.connection_stiffnesses.resize(3);

	mapping.setBeamTopology(topology);

	auto children_of_0 = mapping.getBeamTopology().getChildren(0);
	EXPECT_EQ(children_of_0.size(), 2);
	EXPECT_EQ(children_of_0[0], 1);
	EXPECT_EQ(children_of_0[1], 2);

	auto children_of_1 = mapping.getBeamTopology().getChildren(1);
	EXPECT_TRUE(children_of_1.empty());
}
