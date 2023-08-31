//
// Created by younes on 07/06/2021.
//

#include <Cosserat/config.h>

#include <gtest/gtest.h>
#include <sofa/testing/BaseTest.h>
#include <sofa/defaulttype/config.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/simulation/Node.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/simulation/graph/DAGSimulation.h>

#include <sofa/simulation/graph/SimpleApi.h>
#include <sofa/simulation/common/SceneLoaderXML.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/system/PluginManager.h>

#include <Cosserat/forcefield/BeamHookeLawForceField.inl>
#include <sofa/testing/NumericTest.h>

using sofa::testing::BaseTest ;
using testing::Test;
using sofa::simulation::SceneLoaderXML ;
using namespace sofa::simpleapi;

namespace sofa {

template <typename _DataTypes>
struct BeamHookeLawForceFieldTest : public testing::NumericTest<> {
    typedef _DataTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename Coord::value_type Real;



    typedef sofa::component::forcefield::BeamHookeLawForceField<DataTypes> TheBeamHookeLawForceField;

    // Sets up the test fixture.
    void SetUp() override {
        // initialization or some code to run before each test
        fprintf(stderr, "Starting up ! \n");
        sofa::simpleapi::importPlugin("Sofa.Component");
        sofa::simpleapi::importPlugin("Cosserat");
    }

    // Tears down the test fixture.
    void TearDown() override {
        // code to run after each test;
        // can be used instead of a destructor,
        // but exceptions can be handled in this function only
        if(root) {
            sofa::simulation::node::unload(root);
        }
        fprintf(stderr, "Starting down ! \n");
    }

    void reinit()
    {
        /* @todo do the code here in order to test the updateList function */
        //        m_constraint->UpdateList();
        //        m_constraint.UpdateList();
        //typename CosseratConstraint::SPtr constraint = sofa::core::objectmodel::New<CosseratConstraint>();
        //        simulation::Node::SPtr root = simulation->createNewGraph("root");
        //        root->setGravity( defaulttype::Vector3(0,0,0) );
        auto m_forcefield = new TheBeamHookeLawForceField;

        if (m_forcefield == NULL)
            return ;
        else
            m_forcefield->reinit();
    }

    /**
     *
     */
    void basicAttributesTest();
    void triangle();

//    void addForceTest(const MechanicalParams* mparams,
//                  DataVecDeriv& f ,
//                  const DataVecCoord& x ,
//                  const DataVecDeriv& v) override;
//
//    void addDForceTest(const MechanicalParams* mparams,
//                   DataVecDeriv&   df ,
//                   const DataVecDeriv&
//                   dx ) override;
//
//
//    void addKToMatrixTest(const MechanicalParams* mparams,
//                      const MultiMatrixAccessor* matrix) override;
//
//    double getPotentialEnergyTest(const MechanicalParams* mparams,
//                              const DataVecCoord& x) const override;

protected:
    ///< Root of the scene graph, created by the constructor an re-used in the tests
    simulation::Node::SPtr root;

    void testFonctionnel();
};


template<>
void BeamHookeLawForceFieldTest<defaulttype::Vec3Types>::testFonctionnel() {
    EXPECT_MSG_NOEMIT(Error, Warning) ;
    createObject(root, "MechanicalObject", {{"position", "-1 0 1  1 0 1  -1 0 -1  1 0 -1  0 0 1  0 0 -1  -1 0 0  1 0 0  0 0 0"}});
    createObject(root, "TriangleSetTopologyContainer", {{"triangles", "7 5 8  8 2 6  4 6 0  1 8 4  7 3 5  8 5 2  4 8 6  1 7 8"}});

    auto traction = dynamic_cast<const TheBeamHookeLawForceField *>(
            createObject(root, "BeamHookeLawForceField", {{"name", "beamHookeLaw"},
                                                          {"crossSectionShape", "circular"},
                                                          {"lengthY", "35e-5"},
                                                          {"lengthZ", "1374e-5"},
                                                          {"radius", "0.25"},
                                                          {"varianteSections", "true"}}).get()
    );

    EXPECT_NE(traction, nullptr);
    EXPECT_NE(root.get(), nullptr) ;
    root->init(sofa::core::execparams::defaultInstance()) ;

    auto total_load = dynamic_cast<sofa::core::objectmodel::Data<double> *>(traction->findData("lengthY"));
    for (unsigned int step = 1; step <= 5; ++step) {
        sofa::simulation::node::animate(root.get(), 1);
        EXPECT_DOUBLE_EQ(total_load->getValue(), 4*step) << "Total load at time step " << step << " is incorrect.";
    }
}



template<>
void BeamHookeLawForceFieldTest<defaulttype::Vec3Types>::basicAttributesTest(){
    EXPECT_MSG_NOEMIT(Error) ;

    std::stringstream scene ;
    scene << "<?xml version='1.0'?>"
             "<Node 	name='Root' gravity='0 -9.81 0' time='0' animate='0' >              \n"
             "   <DefaultAnimationLoop/>                                                    \n"
             "   <MechanicalObject name='mstate' template='"<<  DataTypes::Name() << "'/>   \n"
             "   <BeamHookeLawForceField name='myPlaneForceField'/>                         \n"
             "   </Node>                                                                    \n" ;

    Node::SPtr root = SceneLoaderXML::loadFromMemory ("testscene",
                                                      scene.str().c_str()) ;

    EXPECT_NE(root.get(), nullptr) ;
    root->init(sofa::core::execparams::defaultInstance()) ;

    TheBeamHookeLawForceField* forcefield ;
    root->getTreeObject(forcefield) ;

    EXPECT_NE(nullptr, forcefield) ;

    /// List of the supported attributes the user expect to find
    /// This list needs to be updated if you add an attribute.
    sofa::type::vector<std::string> attrnames = {
            "crossSectionShape","youngModulus","poissonRatio","length", "radius",
            "innerRadius", "lengthY", "lengthZ", "varianteSections", "youngModululsList", "poissonRatioList"
    };

    for(auto& attrname : attrnames)
        EXPECT_NE( nullptr, forcefield->findData(attrname) )
        << "Missing attribute with name '"
        << attrname << "'." ;

    for(int i=0; i<10; i++){
        sofa::simulation::node::animate(root.get(),(double)0.01);
    }
}


/***
 *  The test section
 */

// Define the list of DataTypes to instanciate
using ::testing::Types;
typedef Types<defaulttype::Vec3Types> DataTypes; // the types to instantiate.

// Test suite for all the instanciations
TYPED_TEST_SUITE(BeamHookeLawForceFieldTest, DataTypes);// first test case
TYPED_TEST( BeamHookeLawForceFieldTest , basicAttributesTest )
{
    ASSERT_NO_THROW (this->basicAttributesTest());
}

TYPED_TEST(BeamHookeLawForceFieldTest, testFonctionnel) {
        ASSERT_NO_THROW (this->testFonctionnel());
}





}
