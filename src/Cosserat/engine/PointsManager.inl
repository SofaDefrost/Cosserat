#pragma once

#include <Cosserat/config.h>
#include <Cosserat/engine/PointsManager.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseTopology.h>
#include <sofa/component/topology/container/dynamic/PointSetTopologyContainer.h>
#include <sofa/component/topology/container/dynamic/PointSetTopologyModifier.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyContainer.h>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/behavior/Constraint.h>

#include <sofa/type/Quat.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/component/statecontainer/MechanicalObject.h>

namespace cosserat::controller
{
using sofa::core::objectmodel::KeypressedEvent;

PointsManager::PointsManager()
    : d_beamTip(initData(&d_beamTip, "beamTip", "The beam tip")),
      d_radius(initData(&d_radius, double(1), "radius", "sphere radius")),
      d_color(initData(&d_color, sofa::type::Vec4f(1, 0, 0, 1), "color", "Default color is (1,0,0,1)")),
      d_beamPath(initData(&d_beamPath, "beamPath", "path to beam state"))
{
    this->f_listening.setValue(true);
}

PointsManager::~PointsManager(){}

void PointsManager::init()
{
    Inherit1::init();

    if (getTopology() == NULL)
        msg_error() << "Error cannot find the topology";

    if (getMstate() == NULL)
        msg_error() << "Error cannot find the mechanical state";

    this->getContext()->get(m_beam, d_beamPath.getValue());
    if (m_beam == nullptr)
        msg_error() << "Cannot find the beam collision state : " << d_beamPath.getValue();

    this->getContext()->get(m_modifier);
    if (m_modifier == NULL)
    {
        msg_error() << " Error cannot find the EdgeSetTopologyModifier";
        return;
    }
}

void PointsManager::addNewPointToState()
{
    auto x = sofa::helper::getWriteAccessor(*getMstate()->write(sofa::core::VecCoordId::position()));
    auto xRest = sofa::helper::getWriteAccessor(*getMstate()->write(sofa::core::VecCoordId::restPosition()));
    auto xfree = sofa::helper::getWriteAccessor(*getMstate()->write(sofa::core::VecCoordId::freePosition()));
    auto xforce = sofa::helper::getWriteAccessor(*getMstate()->write(sofa::core::VecDerivId::force()));
    const auto &beam = m_beam->readPositions();

    unsigned nbPoints = this->getTopology()->getNbPoints();
    // do not take the last point because there is a bug
    // TODO(dmarchal 2024-07-12) fix the bug you are refering to in the previous line

    size_t beamSz = beam.size();

    m_modifier->addPoints(1, true);
    sofa::type::Vec3 pos = beam[beamSz - 1];

    x.resize(nbPoints + 1);
    xRest.resize(nbPoints + 1);
    xfree.resize(nbPoints + 1);
    xforce.resize(nbPoints + 1);

    x[nbPoints] = pos;
    xRest[nbPoints] = pos;
    xfree[nbPoints] = pos;
    xforce[nbPoints] = sofa::type::Vec3(0, 0, 0);

    m_modifier->notifyEndingEvent();
}

void PointsManager::removeLastPointfromState()
{
    // do not take the last point because there is a bug
    // TODO(dmarchal 2024-07-12) fix the bug you are refering to in the previous line

    unsigned nbPoints = getTopology()->getNbPoints();

    if (nbPoints == 0)
    {
        msg_info() << "No more points";
        return;
    }

    auto x = sofa::helper::getWriteAccessor(*getMstate()->write(sofa::core::VecCoordId::position()));
    auto xfree = sofa::helper::getWriteAccessor(*getMstate()->write(sofa::core::VecCoordId::freePosition()));

    sofa::type::vector<unsigned int> indices = {nbPoints - 1};
    m_modifier->removePoints(indices, true);

    x.resize(nbPoints - 1);
    msg_info() << "the size is equal :" << nbPoints;
    xfree.resize(nbPoints - 1);
    m_modifier->notifyEndingEvent();
}

void PointsManager::handleEvent(sofa::core::objectmodel::Event *event)
{
    if (KeypressedEvent::checkEventType(event))
    {
        KeypressedEvent *ev = static_cast<KeypressedEvent *>(event);
        switch (ev->getKey())
        {
        case 'S':
        case 's':
            msg_info() << "A point is created ."  ;
            addNewPointToState();
            break;
        case 'L':
        case 'l':
            msg_info() <<("Remove point from state \n");
            removeLastPointfromState();
            break;
        }
    }
}

} // namespace cosserat::controller
