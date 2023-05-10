//
// Created by younes on 02/06/2021.
//

#include <gtest/gtest.h>
using testing::Test;

#include <sofa/testing/BaseTest.h>
using sofa::testing::BaseTest ;

#include <SofaTest/Sofa_test.h>
#include <sofa/defaulttype/config.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/simulation/Node.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <SofaSimulationGraph/DAGSimulation.h>

#include <SofaSimulationGraph/SimpleApi.h>
#include <SofaSimulationCommon/SceneLoaderXML.h>
#include <sofa/helper/logging/Message.h>
#include "../../src/constraint/CosseratUnilateralInteractionConstraint.h"
#include "../../src/constraint/CosseratUnilateralInteractionConstraint.inl"

#include "Constraint.h"
#include "../Example.h"

using sofa::simulation::SceneLoaderXML ;

namespace sofa {

template <typename _DataTypes>
struct CosseratUnilateralInteractionConstraintTest : public NumericTest<>
{
    typedef _DataTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename Coord::value_type Real;
    typedef sofa::component::constraintset::CosseratUnilateralInteractionConstraint<DataTypes> CosseratUnilateralInteractionConstraint;

    int* x;
    int GetX() const{
        return *x;
    }

    // Sets up the test fixture.
    void SetUp() override {
        // initialization or some code to run before each test
        fprintf(stderr, "Starting up ! \n");
        x = new int(42);

        sofa::simpleapi::importPlugin("Sofa.Component");
        sofa::simpleapi::importPlugin("Cosserat");

        //create the context for
        if(simulation==nullptr)
            sofa::simulation::setSimulation(simulation = new sofa::simulation::graph::DAGSimulation());
    }

    // Tears down the test fixture.
    void TearDown() override {
        // code to run after each test;
        // can be used instead of a destructor,
        // but exceptions can be handled in this function only
        fprintf(stderr, "Starting down ! \n");
        delete x;
        if(root)
            simulation->unload(root);
    }

    bool UpdateListTest()
    {
        /* @todo do the code here in order to test the updateList function */
        //        m_constraint->UpdateList();
        //        m_constraint.UpdateList();
        //typename CosseratConstraint::SPtr constraint = sofa::core::objectmodel::New<CosseratConstraint>();
        //        simulation::Node::SPtr root = simulation->createNewGraph("root");
        //        root->setGravity( defaulttype::Vector3(0,0,0) );
        CosseratUnilateralInteractionConstraint * m_constraint = new CosseratUnilateralInteractionConstraint;

        if (m_constraint == NULL)
            return false;
        else
            return m_constraint->UpdateList();
    }

    /**
     *
     */
    void attributesTests();

    protected:
    ///< Root of the scene graph, created by the constructor an re-used in the tests
    simulation::Node::SPtr root;
    ///< created by the constructor an re-used in the tests
    simulation::Simulation* simulation {nullptr};
};


template<>
void CosseratUnilateralInteractionConstraintTest<defaulttype::Vec3Types>::attributesTests(){
    /// I'm using '\n' so that the XML parser correctly report the line number
    /// in case of problems.
    std::stringstream scene;
    scene << "<?xml version='1.0'?>                                       \n"
             "<Node name='Root' gravity='0 0 0' time='0' animate='0'   > \n"
             "<MechanicalObject name='o1' template='"<< DataTypes::Name() << "' position='1 2 3'/> \n"
             "<CosseratUnilateralInteractionConstraint template='"<< DataTypes::Name() << "' /> \n"
             "</Node>                                                     \n" ;

    Node::SPtr root = SceneLoaderXML::loadFromMemory ("testscene",
                                                      scene.str().c_str(),
                                                      scene.str().size()) ;
    root->init(sofa::core::execparams::defaultInstance()) ;
    CosseratUnilateralInteractionConstraint * constraint = root->getTreeObject<CosseratUnilateralInteractionConstraint>() ;

        EXPECT_TRUE( constraint != nullptr ) ;
        EXPECT_TRUE( constraint->findData("value") != nullptr ) ;
        EXPECT_TRUE( constraint->findData("force_dumping") != nullptr ) ;
        EXPECT_TRUE( constraint->findData("valueIndex") != nullptr ) ;
        EXPECT_TRUE( constraint->findData("entryPoint") != nullptr ) ;
        EXPECT_TRUE( constraint->findData("vectorOfIndices") != nullptr ) ;

        return ;
    }


/* an easy example */
//TEST(ConstraintTest, Square){
//    int x = 5;
//    int expectedSquare = x*x;
//    EXPECT_EQ(expectedSquare, Square(x));
//}


// Define the list of DataTypes to instantiation
using ::testing::Types;
typedef Types<defaulttype::Vec3Types> DataTypes; // the types to instantiate.

// Test suite for Vec3Types instantiation
TYPED_TEST_SUITE(CosseratUnilateralInteractionConstraintTest, DataTypes);

//TYPED_TEST(ConstraintTest,  MAC)
//{
//    int y = 16;
//    int sum = 100;
//    auto val = this->GetX() ;
////    auto val = ConstraintTest<Vec3Types>->GetX()
//    auto oldSum = sum + val * y; //sum + this->GetX() * y;
//    int expectedNewSum = sum + val * y; //sum + this->GetX() * y;
//
//    EXPECT_EQ(expectedNewSum,MAC(val,y,sum));
//    //EXPECT_EQ(expectedNewSum,oldSum);
//}

// first test case, test the UpdateList function
TYPED_TEST( CosseratUnilateralInteractionConstraintTest, UpdateListTest )
{
    EXPECT_TRUE(this->UpdateListTest());
}

TYPED_TEST( CosseratUnilateralInteractionConstraintTest ,  attributesTests)
{
    ASSERT_NO_THROW(  this->attributesTests() );
}

} //sofa
