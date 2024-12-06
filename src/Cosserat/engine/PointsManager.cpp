
#include "PointsManager.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <math.h>
#include <assert.h> /* assert */

namespace Cosserat
{

    using namespace sofa::defaulttype;

    SOFA_DECL_CLASS(PointsManager)

    void registerPointsManager(sofa::core::ObjectFactory* factory)
    {
        factory->registerObjects(sofa::core::ObjectRegistrationData("add and remove "
                                                        "points from the state will taking into account the changes inside the modifier and the container")
                                      .add<sofa::core::behavior::PointsManager>());
    }

} // namespace sofa
