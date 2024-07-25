#pragma once

#include <Cosserat/config.h>
#include <sofa/core/behavior/BaseController.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/component/topology/container/dynamic/fwd.h>
#include <sofa/core/behavior/MechanicalState.h>

namespace cosserat::controller
{
namespace {
using sofa::component::topology::container::dynamic::PointSetTopologyModifier;
using sofa::core::topology::TopologyContainer;
}

// TODO(dmarchal: 2024-07-12) This class looks like a Controller (not a BaseObject)
class SOFA_COSSERAT_API PointsManager : public sofa::core::behavior::BaseController
{
public:
    SOFA_CLASS(PointsManager, sofa::core::behavior::BaseController);

    PointsManager();
    ~PointsManager();

    sofa::Data<sofa::type::Vec3> d_beamTip;
    sofa::Data<double> d_radius;
    sofa::Data<sofa::type::Vec4f> d_color;
    sofa::Data<std::string> d_beamPath;

    // Inherited from BaseObject
    void init() override;
    void handleEvent(sofa::core::objectmodel::Event *event) override;

    void addNewPointToState();
    void removeLastPointfromState();

private:
    using MState = sofa::core::behavior::MechanicalState<sofa::defaulttype::Vec3Types>;

    PointSetTopologyModifier *m_modifier;
    MState *m_beam;


    auto getTopology() -> TopologyContainer*
    {
        return dynamic_cast<TopologyContainer*>(getContext()->getTopology());
    }

    auto getMstate() -> MState*
    {
        return dynamic_cast<MState *>(getContext()->getMechanicalState());
    }
};

} // namespace cosserat::controller
