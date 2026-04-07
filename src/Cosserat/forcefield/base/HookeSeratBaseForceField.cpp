#include <Cosserat/forcefield/base/HookeSeratBaseForceField.inl>
#include <sofa/defaulttype/VecTypes.h>
namespace sofa::component::forcefield {

	template class SOFA_COSSERAT_API HookeSeratBaseForceField<defaulttype::Vec3Types>;
	template class SOFA_COSSERAT_API HookeSeratBaseForceField<defaulttype::Vec6Types>;

} // namespace sofa::component::forcefield
