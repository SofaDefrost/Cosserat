#ifndef DrawTrianglesComponent_H
#define DrawTrianglesComponent_H

#include <string>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/helper/system/config.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/visual/DrawTool.h>


#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/helper/ColorMap.h>
#include <SofaSimpleFem/TetrahedronFEMForceField.h>


namespace sofa {

namespace component {

namespace controller {

using namespace defaulttype;

class DrawTrianglesComponent :  public core::objectmodel::BaseObject
{
public:

    SOFA_CLASS(DrawTrianglesComponent,core::objectmodel::BaseObject);
    typedef typename defaulttype::Vec3dTypes DataTypes;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::VecReal VecReal;
    typedef VecCoord Vector;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;

    typedef core::objectmodel::Data<VecDeriv>    DataVecDeriv;
    typedef core::objectmodel::Data<VecCoord>    DataVecCoord;

    typedef core::topology::BaseMeshTopology::index_type Index;
    typedef core::topology::BaseMeshTopology::Tetra Element;
    typedef core::topology::BaseMeshTopology::SeqTetrahedra VecElement;
    typedef core::topology::BaseMeshTopology::Tetrahedron Tetrahedron;
    typedef core::topology::BaseMeshTopology::Triangle Triangle;

    DrawTrianglesComponent();

    void init();
    void draw(const core::visual::VisualParams* vparams);
    void handleEvent(sofa::core::objectmodel::Event* event);

    Data<double>                               d_transparency;
    Data<helper::vector<Triangle> >            d_vecTriangles;
    Data<helper::vector<Tetrahedron> >         d_vecTetra;
    Data<double>                               d_maxStress;
    Data<double>                               d_minStress;
    Data<helper::vector<double>>               d_maxVonMisesPerNode;
    Data<bool>                                 d_draw;

    helper::ColorMap m_VonMisesColorMap;
    sofa::component::forcefield::TetrahedronFEMForceField<Vec3dTypes>* m_tetraForceField;
    core::behavior::MechanicalState<Vec3dTypes>* m_state;

    double m_minVM, m_maxVM;

};

} //end namespace controller

} //end namespace component

} //end namespace sofa

#endif // DrawTrianglesComponent_H
