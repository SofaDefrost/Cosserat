#pragma once

#include <Cosserat/config.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/behavior/BaseController.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/component/topology/container/dynamic/PointSetTopologyContainer.h>
#include <sofa/component/topology/container/dynamic/PointSetTopologyModifier.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyContainer.h>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>
#include <sofa/helper/AdvancedTimer.h>
// #include <sofa/gl/template.h>
#include <sofa/core/behavior/Constraint.h>

typedef sofa::component::topology::container::dynamic::PointSetTopologyModifier PointSetTopologyModifier;

using sofa::core::objectmodel::KeypressedEvent;
namespace sofa::core::behavior
{
    class SOFA_COSSERAT_API PointsManager : public sofa::core::objectmodel::BaseObject
    {

    public:
        SOFA_CLASS(PointsManager, sofa::core::objectmodel::BaseObject);

        typedef sofa::defaulttype::Vec3dTypes DataTypes;
        typedef DataTypes::VecCoord VecCoord;
        typedef DataTypes::Coord Coord;
        typedef DataTypes::Real Real;

    public:
        PointsManager();
        typedef type::Vec3 Vec3;

        Data<Vec3> d_beamTip;
        Data<double> d_radius;
        Data<type::Vec4f> d_color;
        Data<std::string> d_beamPath;

        PointSetTopologyModifier *m_modifier;
        core::behavior::MechanicalState<DataTypes> *m_beam;

        void init() override;
        void handleEvent(sofa::core::objectmodel::Event *event) override;
        // void draw(const core::visual::VisualParams *vparams);

        void addNewPointToState();
        void removeLastPointfromState();

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
