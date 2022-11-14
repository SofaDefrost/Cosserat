#pragma once

#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseTopology/PointSetTopologyModifier.h>
#include <SofaBaseTopology/PointSetTopologyModifier.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/behavior/BaseController.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseTopology.h>
#include <SofaBaseTopology/PointSetTopologyContainer.h>
#include <SofaBaseTopology/EdgeSetTopologyContainer.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <sofa/helper/AdvancedTimer.h>
// #include <sofa/gl/template.h>
#include <sofa/core/behavior/Constraint.h>

using sofa::core::objectmodel::KeypressedEvent;
namespace sofa::core::behavior
{
    class PointsManager : public sofa::core::objectmodel::BaseObject
    {

    public:
        SOFA_CLASS(PointsManager, sofa::core::objectmodel::BaseObject);

        typedef sofa::defaulttype::Vec3dTypes DataTypes;
        typedef DataTypes::VecCoord VecCoord;
        typedef DataTypes::Coord Coord;
        typedef DataTypes::Real Real;

    public:
        PointsManager();
        typedef type::Vector3 Vector3;

        Data<Vector3> d_beamTip;
        Data<double> d_radius;
        Data<type::Vec4f> d_color;
        Data<std::string> d_beamPath;

        sofa::component::topology::PointSetTopologyModifier *m_modifier;
        core::behavior::MechanicalState<DataTypes> *m_beam;

        void init();
        void handleEvent(sofa::core::objectmodel::Event *event);
        // void draw(const core::visual::VisualParams *vparams);

        virtual void addNewPointToState();
        virtual void removeLastePointfromState();

        topology::TopologyContainer *getTopology()
        {
            return dynamic_cast<topology::TopologyContainer *>(getContext()->getTopology());
        }

        sofa::core::behavior::MechanicalState<DataTypes> *getMstate()
        {
            return dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes> *>(getContext()->getMechanicalState());
        }
    };

} // namespace sofa
