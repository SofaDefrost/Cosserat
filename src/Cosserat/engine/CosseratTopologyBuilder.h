#pragma once

#include <Cosserat/engine/CosseratBeamGeometryEngine.h>
#include <Cosserat/engine/CosseratIntrinsicState.h>
#include <liegroups/SO3.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyContainer.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/SingleLink.h>
#include <sofa/type/Vec.h>

namespace sofa::component::cosserat::engine {

using Vec3dTypes = sofa::defaulttype::Vec3dTypes;
using SO3 = sofa::component::cosserat::liegroups::SO3<double>;

// Note: Kept for compatibility if used elsewhere, but not needed for IntrinsicState
struct SO3Types {
    using Coord = SO3;
    using Deriv = SO3::TangentVector;
    using VecCoord = sofa::helper::vector<Coord>;
    using VecDeriv = sofa::helper::vector<Deriv>;
    using Real = double;
};

/**
 * @brief Builder for initializing a Cosserat Topology and States
 */
class CosseratTopologyBuilder : public sofa::core::objectmodel::BaseObject {
   public:
    SOFA_CLASS(CosseratTopologyBuilder, sofa::core::objectmodel::BaseObject);

    CosseratTopologyBuilder();
    virtual ~CosseratTopologyBuilder() = default;

    void init() override;

    // Data for physical parameters
    sofa::core::objectmodel::Data<double> d_totalLength;
    sofa::core::objectmodel::Data<int> d_nbSections;
    sofa::core::objectmodel::Data<double> d_radius;
    sofa::core::objectmodel::Data<double> d_youngModulus;
    sofa::core::objectmodel::Data<double> d_shearModulus;

    // Data for rest states (stiffness)
    sofa::core::objectmodel::Data<double> d_EA;
    sofa::core::objectmodel::Data<double> d_GA;
    sofa::core::objectmodel::Data<double> d_EIy;
    sofa::core::objectmodel::Data<double> d_EIz;
    sofa::core::objectmodel::Data<double> d_GJ;

    // Links to target components
    sofa::core::objectmodel::SingleLink<CosseratIntrinsicState> l_intrinsicState;
    sofa::core::objectmodel::SingleLink<
        sofa::component::topology::container::dynamic::EdgeSetTopologyContainer>
        l_topology;
};

}  // namespace sofa::component::cosserat::engine
