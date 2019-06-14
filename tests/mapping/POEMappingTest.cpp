#include <string>
using std::string ;
#include <SofaTest/Multi2Mapping_test.h>
#include <SofaSimulationGraph/DAGSimulation.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/simulation/UpdateLinksVisitor.h>
#include <sofa/simulation/InitVisitor.h>
#include <sstream>
#include <sofa/core/MechanicalParams.h>
#include <sofa/simulation/VectorOperations.h>
#include <SofaBaseLinearSolver/FullVector.h>
#include <SofaEigen2Solver/EigenSparseMatrix.h>
#include <SofaBaseMechanics/MechanicalObject.h>

#include "../../src/mapping/POEMapping.inl"
#include <SceneCreator/SceneCreator.h>
#include <sofa/helper/vector.h>
#include <sofa/core/MultiMapping.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/simulation/SceneLoaderFactory.h>

#include <SofaSimulationCommon/SceneLoaderXML.h>
#include <SofaPython/SceneLoaderPY.h>

using sofa::simulation::SceneLoaderXML ;
using sofa::simulation::Node ;
using sofa::simulation::SceneLoaderPY;
using sofa::component::container::MechanicalObject ;

namespace sofa {
namespace { // anonymous namespace
using namespace core;
using namespace component;
using defaulttype::Vec;
using defaulttype::Mat;


/**  Test suite for RigidMapping.
The test cases are defined in the #Test_Cases member group.
  */
template <typename _POEMapping>
struct POEMappingTest : public Multi2Mapping_test<_POEMapping>
{

    typedef _POEMapping POEMapping;
    typedef Multi2Mapping_test<POEMapping> Inherit;

    typedef typename POEMapping::In1 In1Type;
    typedef typename In1Type::VecCoord In1VecCoord;
    typedef typename In1Type::VecDeriv In1VecDeriv;
    typedef typename In1Type::Coord In1Coord;
    typedef typename In1Type::Deriv In1Deriv;
    typedef typename POEMapping::In1DataVecCoord In1DataVecCoord;
    typedef typename POEMapping::In1DataVecDeriv In1DataVecDeriv;
    typedef container::MechanicalObject<In1Type> In1DOFs;
    typedef typename In1DOFs::ReadVecCoord  ReadIn1VecCoord;
    typedef typename In1DOFs::WriteVecCoord WriteIn1VecCoord;
    typedef typename In1DOFs::WriteVecDeriv WriteIn1VecDeriv;
    typedef typename In1Type::Real Real;
    typedef Mat<In1Type::spatial_dimensions, In1Type::spatial_dimensions, Real> In1RotationMatrix;

    typedef typename POEMapping::In2 In2Type;
    typedef typename In2Type::VecCoord In2VecCoord;
    typedef typename In2Type::VecDeriv In2VecDeriv;
    typedef typename In2Type::Coord In2Coord;
    typedef typename In2Type::Deriv In2Deriv;
    typedef typename POEMapping::In2DataVecCoord In2DataVecCoord;
    typedef typename POEMapping::In2DataVecDeriv In2DataVecDeriv;
    typedef container::MechanicalObject<In2Type> In2DOFs;
    typedef typename In2DOFs::ReadVecCoord  ReadIn2VecCoord;
    typedef typename In2DOFs::WriteVecCoord WriteIn2VecCoord;
    typedef typename In2DOFs::WriteVecDeriv WriteIn2VecDeriv;
    typedef typename In1Type::Real Real2;
    typedef Mat<In2Type::spatial_dimensions, In2Type::spatial_dimensions, Real> In2RotationMatrix;

    typedef typename POEMapping::Out OutType;
    typedef typename OutType::VecCoord OutVecCoord;
    typedef typename OutType::VecDeriv OutVecDeriv;
    typedef typename OutType::Coord OutCoord;
    typedef typename OutType::Deriv OutDeriv;
    typedef typename POEMapping::OutDataVecCoord OutDataVecCoord;
    typedef typename POEMapping::OutDataVecDeriv OutDataVecDeriv;
    typedef container::MechanicalObject<OutType> OutDOFs;
    typedef typename OutDOFs::WriteVecCoord WriteOutVecCoord;
    typedef typename OutDOFs::WriteVecDeriv WriteOutVecDeriv;
    typedef typename OutDOFs::ReadVecCoord ReadOutVecCoord;
    typedef typename OutDOFs::ReadVecDeriv ReadOutVecDeriv;

    typedef core::Multi2Mapping <In1Type, In2Type, OutType> Multi2Mapping;

    typedef component::linearsolver::EigenSparseMatrix<In1Type, OutType> SparseJMatrixEigen1;
    typedef component::linearsolver::EigenSparseMatrix<In2Type, OutType> SparseJMatrixEigen2;
    typedef linearsolver::EigenSparseMatrix<In1Type, In2Type> SparseKMatrixEigen1;
    typedef linearsolver::EigenSparseMatrix<In2Type, In2Type> SparseKMatrixEigen2;



    POEMapping* poeMapping; ///< the mapping to be tested
    vector<In1DOFs*>  in1Dofs; ///< mapping input
    vector<In2DOFs*>  in2Dofs; ///< mapping input
    OutDOFs* outDofs; ///< mapping output

    //simulation::SceneLoaderPY loader;

    /// Constructor
    POEMappingTest()
    {
        //        this->errorFactorDJ = 200;

        //        poeMapping = static_cast<POEMapping*>( this->Multi2Mapping );

        //        // beamLengthMapping::getJs is not yet implemented
        //        this->flags &= ~Inherit::TEST_getJs;

        //        // beamLengthMapping::getK is not yet implemented
        //        //this->flags &= ~Inherit::TEST_getK;

        //        // beamLengthMapping::applyDJT is not yet implemented
        //        //this->flags &= ~Inherit::TEST_applyDJT;

    }

    /** @name Test_Cases
      verify the computation of the beam lengths and the derivatives
      */
    ///@{
    /** Two frames are placed + line topology + and the mapping is constructed
    */

    bool runTest1()
    {
        printf("================================> Inside the POEMapping test scene \n");

        const char *filename ;

        SceneLoaderPY loader;
        Node::SPtr root = loader.doLoad(filename);

        //std::cout<<"*******************  Get the mapping ";
        POEMapping* poeMapping;
        this->root->getTreeObject(poeMapping);
        this->mapping = poeMapping;
        this->deltaRange.first= 1;
        this->deltaRange.second= 1000;


        MechanicalObject<defaulttype::Vec3Types>* FromModel1 = nullptr;
        this->root->getTreeObject(FromModel1);
        //this->in1Dofs = FromModel1;
        const Data<In1VecCoord>& dataIn1X = *FromModel1->read(VecCoordId::position());
        const In1VecCoord& x1in = dataIn1X.getValue();

        helper::vector<In1VecCoord> vecIn1;
        vecIn1.push_back(x1in);

        MechanicalObject<defaulttype::Rigid3dTypes>* FromModel2 = nullptr;
        this->root->getTreeObject(FromModel2);
        //this->in2Dofs = FromModel2;
        const Data<In2VecCoord>& dataIn2X = *FromModel2->read(VecCoordId::position());
        const In2VecCoord& x2in = dataIn2X.getValue();

        helper::vector<In2VecCoord> vecIn2;
        vecIn2.push_back(x2in);


        MechanicalObject<defaulttype::Rigid3dTypes>* ToModel = nullptr;
        this->root->getTreeObject(ToModel);
        this->outDofs= ToModel;
        const Data<OutVecCoord>& dataOutX = *ToModel->read(VecCoordId::position());
        const OutVecCoord& xout = dataOutX.getValue();

        std::cout<<" x1 in  = "<<x1in<<std::endl;

        std::cout<<" x2 in  = "<<x2in<<std::endl;

        std::cout<<" x out  = "<<xout<<std::endl;



        // new values of DOFs: exact same position
        //        In1VecCoord x1in_new(Nin);
        //        x1in_new[1][0]=1.0;
        //        x1in_new[2][0]=2.0;

        //        // new values of DOFs: exact same position
        //        In2VecCoord x2in_new(Nin);
        //        x2in_new[1][0]=1.0;
        //        x2in_new[2][0]=2.0;


        OutVecCoord expectedoutcoord;
        OutCoord r0 = OutCoord(defaulttype::Vec3(0.5,0,0),defaulttype::Quat(0.0374912, 0, 0, 0.999297));
        OutCoord r1 = OutCoord(defaulttype::Vec3(0.5,0,0),defaulttype::Quat(0.0374912, 0, 0, 0.999297));
        OutCoord r2 = OutCoord(defaulttype::Vec3(0.5,0,0),defaulttype::Quat(0.0374912, 0, 0, 0.999297));
        OutCoord r3 = OutCoord(defaulttype::Vec3(0.5,0,0),defaulttype::Quat(0.0374912, 0, 0, 0.999297));
        OutCoord r4 = OutCoord(defaulttype::Vec3(0.5,0,0),defaulttype::Quat(0.0374912, 0, 0, 0.999297));
        expectedoutcoord.push_back(r0);
        expectedoutcoord.push_back(r1);
        expectedoutcoord.push_back(r2);
        expectedoutcoord.push_back(r3);
        expectedoutcoord.push_back(r4);

        ///Parent function
        //bool runTest(const vector<In1VecCoord>& in1Coords, const vector<In2VecCoord>& in2Coords, const OutVecCoord& expectedChildCoords)

        //return this->runTest(vecIn1,vecIn2,expectedoutcoord) ;
        return Inherit::runTest(vecIn1,vecIn2,expectedoutcoord);
    }

};

// Define the list of types to instanciate. We do not necessarily need to test all combinations.
using testing::Types;
typedef Types<
mapping::POEMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >> DataTypes; // the types to instanciate.


// Test suite for all the instanciations
TYPED_TEST_CASE(POEMappingTest, DataTypes);
// first test case
TYPED_TEST( POEMappingTest , runTest1 )
{
    ASSERT_TRUE(this->runTest1());
}


}// anonymous namespace
}//sofa
