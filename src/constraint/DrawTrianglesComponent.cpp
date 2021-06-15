#include "DrawTrianglesComponent.inl"
//#include<sofa/helper/system/config.h>

namespace sofa {

namespace component {

namespace controller {

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(DrawTrianglesComponent)

int DrawTrianglesComponentClass = core::RegisterObject("Connection to Motive")
.add< DrawTrianglesComponent >();


} //end namespace controller

} //end namespace component

} //end namespace sofa
