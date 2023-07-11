
#include "PointsManager.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/visual/VisualParams.h>
#include <math.h>
#include <assert.h> /* assert */

namespace sofa::core::behavior
{

    using namespace sofa::defaulttype;

    SOFA_DECL_CLASS(PointsManager)

    int PointsManagerClass = core::RegisterObject("add and remove "
                                                  "points from the state will taking into account the changes inside the modifier and the container")
                                .add<PointsManager>();

} // namespace sofa
