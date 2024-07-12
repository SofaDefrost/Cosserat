
#include <Cosserat/engine/PointsManager.inl>
#include <sofa/core/ObjectFactory.h>

namespace cosserat::controller
{
int PointsManagerClass = sofa::core::RegisterObject(
                             R"doc(Add or remove controls point in a beam structure)doc")
                         .add<PointsManager>();
} // namespace cosserat::controller
