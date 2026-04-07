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

#include <sofa/simpleapi/SimpleApi.h>
#include <sofa/simulation/common/SceneLoaderXML.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/system/PluginManager.h>
#include <sofa/helper/system/FileSystem.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
#include <sofa/simulation/TaskScheduler.h>
#include <chrono>

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

    BeamHookeLawForceFieldTest()
    {
        root = sofa::simpleapi::createRootNode(sofa::simulation::getSimulation(), "root");

        createObject(root, "DefaultAnimationLoop");
        createObject(root, "DefaultVisualManagerLoop");
        sofa::simpleapi::importPlugin("Sofa.Component");
        sofa::simpleapi::importPlugin("Cosserat");
    }



    typedef sofa::component::forcefield::BeamHookeLawForceField<DataTypes> TheBeamHookeLawForceField;

    // Sets up the test fixture.
    void doSetUp() override {
        // initialization or some code to run before each test
        fprintf(stderr, "Starting up ! \n");
    }

    // Tears down the test fixture.
    void doTearDown() override {
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
     *
     */
    void basicAttributesTest();
    void testUniformSection();
    void testVariantSection();
    void testParallelPerformance();
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

    // Helper function to create a basic beam model
    TheBeamHookeLawForceField* createBeamModel(bool variantSections, bool useMultiThreading);
    
    // Compare performance between single-threaded and multi-threaded force computation
    void comparePerformance(int numElements, bool variantSections);

    void testFonctionnel();
};


template<typename DataTypes>
typename BeamHookeLawForceFieldTest<DataTypes>::TheBeamHookeLawForceField* 
BeamHookeLawForceFieldTest<DataTypes>::createBeamModel(bool variantSections, bool useMultiThreading) {
    EXPECT_MSG_NOEMIT(Error, Warning);
    ASSERT_NE(root, nullptr);
    
    // Create mechanical object with positions for testing
    sofa::component::statecontainer::MechanicalObject<DataTypes>* mstate = dynamic_cast<sofa::component::statecontainer::MechanicalObject<DataTypes>*>(
        createObject(root, "MechanicalObject", {
            {"name", "mstate"},
            {"position", "0 0 0  0.1 0 0  0.2 0 0  0.3 0 0  0.4 0 0  0.5 0 0  0.6 0 0  0.7 0 0  0.8 0 0  0.9 0 0"}
        }).get()
    );
    
    // Set rest position to match initial position
    mstate->writeRestPositions(mstate->x);
    
    // Create beam force field
    std::vector<std::pair<std::string, std::string>> attributes = {
        {"name", "beamForceField"},
        {"crossSectionShape", "circular"},
        {"youngModulus", "1e5"},
        {"poissonRatio", "0.3"},
        {"radius", "0.01"},
        {"variantSections", variantSections ? "true" : "false"},
        {"useMultiThreading", useMultiThreading ? "true" : "false"},
        {"length", "0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1"}
    };
    
    // Add variant section data if needed
    if (variantSections) {
        attributes.push_back({"youngModulusList", "1e5 1.1e5 1.2e5 1.3e5 1.4e5 1.5e5 1.6e5 1.7e5 1.8e5"});
        attributes.push_back({"poissonRatioList", "0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38"});
    }
    
    auto forcefield = dynamic_cast<TheBeamHookeLawForceField*>(
        createObject(root, "BeamHookeLawForceField", attributes).get()
    );
    
    EXPECT_NE(forcefield, nullptr);
    root->init(sofa::core::execparams::defaultInstance());
    
    return forcefield;
}

template<typename DataTypes>
void BeamHookeLawForceFieldTest<DataTypes>::comparePerformance(int numElements, bool variantSections) {
    // Reset the root node for this test
    if (root) {
        sofa::simulation::node::unload(root);
    }
    root = sofa::simpleapi::createRootNode(sofa::simulation::getSimulation(), "root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");
    
    // Create a large mechanical object for performance testing
    std::stringstream positionStr;
    std::stringstream lengthStr;
    std::stringstream youngModuliStr;
    std::stringstream poissonRatioStr;
    
    for (int i = 0; i < numElements; i++) {
        positionStr << i/10.0 << " 0 0  ";
        if (i < numElements - 1) {
            lengthStr << "0.1 ";
            youngModuliStr << (1.0 + 0.1*i) << "e5 ";
            poissonRatioStr << (0.3 + 0.01*i) << " ";
        }
    }
    
    // Create mechanical object
    sofa::component::statecontainer::MechanicalObject<DataTypes>* mstate = dynamic_cast<sofa::component::statecontainer::MechanicalObject<DataTypes>*>(
        createObject(root, "MechanicalObject", {
            {"name", "mstate"},
            {"position", positionStr.str()}
        }).get()
    );
    
    // Set rest positions
    mstate->writeRestPositions(mstate->x);
    
    // Apply a small deformation to test forces
    auto x = mstate->x.beginEdit();
    for (size_t i = 0; i < x->size(); i++) {
        (*x)[i][0] += 0.001 * i;  // Apply small deformation
    }
    mstate->x.endEdit();
    
    // Run with single-threading
    auto singleThreadForceField = dynamic_cast<TheBeamHookeLawForceField*>(
        createObject(root, "BeamHookeLawForceField", {
            {"name", "singleThreadFF"},
            {"crossSectionShape", "circular"},
            {"youngModulus", "1e5"},
            {"poissonRatio", "0.3"},
            {"radius", "0.01"},
            {"variantSections", variantSections ? "true" : "false"},
            {"useMultiThreading", "false"},
            {"length", lengthStr.str()}
        }).get()
    );
    
    if (variantSections) {
        singleThreadForceField->findData("youngModulusList")->read(youngModuliStr.str());
        singleThreadForceField->findData("poissonRatioList")->read(poissonRatioStr.str());
    }
    
    root->init(sofa::core::execparams::defaultInstance());
    
    // Measure single-threaded performance
    auto start_single = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 10; i++) {  // Run multiple iterations for more accurate timing
        sofa::simulation::node::animate(root.get(), 0.01);
    }
    auto end_single = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> single_thread_time = end_single - start_single;
    
    // Remove single threaded force field
    root->removeObject(singleThreadForceField);
    
    // Run with multi-threading
    auto multiThreadForceField = dynamic_cast<TheBeamHookeLawForceField*>(
        createObject(root, "BeamHookeLawForceField", {
            {"name", "multiThreadFF"},
            {"crossSectionShape", "circular"},
            {"youngModulus", "1e5"},
            {"poissonRatio", "0.3"},
            {"radius", "0.01"},
            {"variantSections", variantSections ? "true" : "false"},
            {"useMultiThreading", "true"},
            {"length", lengthStr.str()}
        }).get()
    );
    
    if (variantSections) {
        multiThreadForceField->findData("youngModulusList")->read(youngModuliStr.str());
        multiThreadForceField->findData("poissonRatioList")->read(poissonRatioStr.str());
    }
    
    root->init(sofa::core::execparams::defaultInstance());
    
    // Measure multi-threaded performance
    auto start_multi = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 10; i++) {  // Run multiple iterations for more accurate timing
        sofa::simulation::node::animate(root.get(), 0.01);
    }
    auto end_multi = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> multi_thread_time = end_multi - start_multi;
    
    // Calculate speedup
    double speedup = single_thread_time.count() / multi_thread_time.count();
    
    // Output performance results
    std::cout << "Performance comparison for " << (variantSections ? "variant" : "uniform") 
              << " sections with " << numElements << " elements:" << std::endl;
    std::cout << "  Single-threaded time: " << single_thread_time.count() << " ms" << std::endl;
    std::cout << "  Multi-threaded time:  " << multi_thread_time.count() << " ms" << std::endl;
    std::cout << "  Speedup factor:       " << speedup << "x" << std::endl;
    
    // We expect some speedup with multithreading, though the exact amount depends on hardware
    // For a reasonable number of elements, we should see at least some improvement
    EXPECT_GT(speedup, 1.0) << "Multithreading should provide some speedup";
}

template<typename DataTypes>
void BeamHookeLawForceFieldTest<DataTypes>::testUniformSection() {
    EXPECT_MSG_NOEMIT(Error, Warning);
    
    // Create models for testing - one single-threaded and one multi-threaded
    if (root) {
        sofa::simulation::node::unload(root);
    }
    root = sofa::simpleapi::createRootNode(sofa::simulation::getSimulation(), "root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");
    
    // Create a force field with uniform sections using single threading
    auto singleThreadFF = createBeamModel(false, false);
    ASSERT_NE(singleThreadFF, nullptr);
    
    // Create a mechanical object and apply a known deformation
    auto mstate = dynamic_cast<sofa::component::statecontainer::MechanicalObject<DataTypes>*>(
        root->getTreeObject("mstate"));
    ASSERT_NE(mstate, nullptr);
    
    // Apply deformation
    auto x = mstate->x.beginEdit();
    (*x)[4][0] += 0.01;  // Apply deformation to middle node
    mstate->x.endEdit();
    
    // Compute forces with single-threading
    sofa::simulation::node::animate(root.get(), 0.01);
    
    // Store the forces from single-threaded computation
    const auto& singleThreadForces = mstate->f.getValue();
    
    // Remove single-threaded force field
    root->removeObject(singleThreadFF);
    
    // Create a force field with uniform sections using multi-threading
    auto multiThreadFF = createBeamModel(false, true);
    ASSERT_NE(multiThreadFF, nullptr);
    
    // Compute forces with multi-threading
    sofa::simulation::node::animate(root.get(), 0.01);
    
    // Compare forces from single-threaded and multi-threaded computations
    const auto& multiThreadForces = mstate->f.getValue();
    
    // The forces should be identical regardless of threading
    ASSERT_EQ(singleThreadForces.size(), multiThreadForces.size());
    for (size_t i = 0; i < singleThreadForces.size(); i++) {
        for (size_t j = 0; j < 3; j++) {  // Compare each component (assuming Vec3)
            EXPECT_NEAR(singleThreadForces[i][j], multiThreadForces[i][j], 1e-10) 
                << "Force difference at index " << i << ", component " << j;
        }
    }
}

template<typename DataTypes>
void BeamHookeLawForceFieldTest<DataTypes>::testVariantSection() {
    EXPECT_MSG_NOEMIT(Error, Warning);
    
    // Create models for testing - one single-threaded and one multi-threaded
    if (root) {
        sofa::simulation::node::unload(root);
    }
    root = sofa::simpleapi::createRootNode(sofa::simulation::getSimulation(), "root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");
    
    // Create a force field with variant sections using single threading
    auto singleThreadFF = createBeamModel(true, false);
    ASSERT_NE(singleThreadFF, nullptr);
    
    // Create a mechanical object and apply a known deformation
    auto mstate = dynamic_cast<sofa::component::statecontainer::MechanicalObject<DataTypes>*>(
        root->getTreeObject("mstate"));
    ASSERT_NE(mstate, nullptr);
    
    // Apply deformation
    auto x = mstate->x.beginEdit();
    for (size_t i = 0; i < x->size(); i++) {
        (*x)[i][0] += 0.001 * i;  // Apply varying deformation
    }
    mstate->x.endEdit();
    
    // Compute forces with single-threading
    sofa::simulation::node::animate(root.get(), 0.01);
    
    // Store the forces from single-threaded computation
    const auto& singleThreadForces = mstate->f.getValue();
    
    // Remove single-threaded force field
    root->removeObject(singleThreadFF);
    
    // Create a force field with variant sections using multi-threading
    auto multiThreadFF = createBeamModel(true, true);
    ASSERT_NE(multiThreadFF, nullptr);
    
    // Compute forces with multi-threading
    sofa::simulation::node::animate(root.get(), 0.01);
    
    // Compare forces from single-threaded and multi-threaded computations
    const auto& multiThreadForces = mstate->f.getValue();
    
    // The forces should be identical regardless of threading
    ASSERT_EQ(singleThreadForces.size(), multiThreadForces.size());
    for (size_t i = 0; i < singleThreadForces.size(); i++) {
        for (size_t j = 0; j < 3; j++) {  // Compare each component (assuming Vec3)
            EXPECT_NEAR(singleThreadForces[i][j], multiThreadForces[i][j], 1e-10) 
                << "Force difference at index " << i << ", component " << j;
        }
    }
}

template<typename DataTypes>
void BeamHookeLawForceFieldTest<DataTypes>::testParallelPerformance() {
    EXPECT_MSG_NOEMIT(Error, Warning);
    
    // Test performance with different numbers of elements
    comparePerformance(100, false);  // Uniform sections with 100 elements
    comparePerformance(500, false);  // Uniform sections with 500 elements
    comparePerformance(100, true);   // Variant sections with 100 elements
    comparePerformance(500, true);   // Variant sections with 500 elements
}

template<typename DataTypes>
void BeamHookeLawForceFieldTest<DataTypes>::testFonctionnel() {
    EXPECT_MSG_NOEMIT(Error, Warning);
    ASSERT_NE(root, nullptr);
    
    // Reset root for this test
    if (root) {
        sofa::simulation::node::unload(root);
    }
    root = sofa::simpleapi::createRootNode(sofa::simulation::getSimulation(), "root");
    createObject(root, "DefaultAnimationLoop");
    createObject(root, "DefaultVisualManagerLoop");
    
    // Create a mechanical object with positions for beam simulation
    sofa::component::statecontainer::MechanicalObject<DataTypes>* mstate = dynamic_cast<sofa::component::statecontainer::MechanicalObject<DataTypes>*>(
        createObject(root, "MechanicalObject", {
            {"name", "mstate"},
            {"position", "0 0 0  0.1 0 0  0.2 0 0  0.3 0 0  0.4 0 0"}
        }).get()
    );
    
    // Set rest position to match initial position
    mstate->writeRestPositions(mstate->x);
    
    // Create beam force field with both uniform and variant options for testing
    auto forcefield = dynamic_cast<TheBeamHookeLawForceField*>(
        createObject(root, "BeamHookeLawForceField", {
            {"name", "beamForceField"},
            {"crossSectionShape", "circular"},
            {"youngModulus", "1e5"},
            {"poissonRatio", "0.3"},
            {"radius", "0.01"},
            {"variantSections", "false"},
            {"useMultiThreading", "true"},
            {"length", "0.1 0.1 0.1 0.1"}
        }).get()
    );
    
    EXPECT_NE(forcefield, nullptr);
    root->init(sofa::core::execparams::defaultInstance());
    
    // Apply deformation and verify forces
    auto x = mstate->x.beginEdit();
    (*x)[2][0] += 0.01;  // Apply deformation to middle node
    mstate->x.endEdit();
    
    // Run one step of simulation
    sofa::simulation::node::animate(root.get(), 0.01);
    
    // Check that forces are computed and are non-zero
    const auto& forces = mstate->f.getValue();
    ASSERT_EQ(forces.size(), 5);
    
    // Verify that forces are computed correctly (non-zero at deformed points)
    bool foundNonZeroForce = false;
    for (size_t i = 0; i < forces.size(); i++) {
        if (forces[i].norm() > 1e-10) {
            foundNonZeroForce = true;
            break;
        }
    }
    
    EXPECT_TRUE(foundNonZeroForce) << "Expected non-zero forces due to deformation";
    
    // Now switch to variant sections and test again
    root->removeObject(forcefield);
    
    auto variantForcefield = dynamic_cast<TheBeamHookeLawForceField*>(
        createObject(root, "BeamHookeLawForceField", {
            {"name", "beamForceField"},
            {"crossSectionShape", "circular"},
            {"youngModulus", "1e5"},
            {"poissonRatio", "0.3"},
            {"radius", "0.01"},
            {"variantSections", "true"},
            {"youngModulusList", "1e5 1.1e5 1.2e5 1.3e5"},
            {"poissonRatioList", "0.3 0.31 0.32 0.33"},
            {"useMultiThreading", "true"},
            {"length", "0.1 0.1 0.1 0.1"}
        }).get()
    );
    
    EXPECT_NE(variantForcefield, nullptr);
    root->init(sofa::core::execparams::defaultInstance());
    
    // Run simulation with variant sections
    sofa::simulation::node::animate(root.get(), 0.01);
    
    // Check that forces are computed and are non-zero with variant sections
    const auto& variantForces = mstate->f.getValue();
    ASSERT_EQ(variantForces.size(), 5);
    
    // Verify that forces are computed correctly with variant sections
    foundNonZeroForce = false;
    for (size_t i = 0; i < variantForces.size(); i++) {
        if (variantForces[i].norm() > 1e-10) {
            foundNonZeroForce = true;
            break;
        }
    }
    
    EXPECT_TRUE(foundNonZeroForce) << "Expected non-zero forces with variant sections";
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
            "innerRadius", "lengthY", "lengthZ", "variantSections", "youngModulusList", "poissonRatioList"
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

        ASSERT_NO_THROW (this->testFonctionnel());
}





}
