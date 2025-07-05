#include <Cosserat/forcefield/standard/BaseBeamForceField.inl>
#include <sofa/defaulttype/VecTypes.h>
namespace sofa::component::forcefield {

	template class SOFA_COSSERAT_API BaseBeamForceField<defaulttype::Vec3Types>;
	template class SOFA_COSSERAT_API BaseBeamForceField<defaulttype::Vec6Types>;

} // namespace sofa::component::forcefield

