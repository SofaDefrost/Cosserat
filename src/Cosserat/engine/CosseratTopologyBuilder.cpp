#include <Cosserat/engine/CosseratTopologyBuilder.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace Cosserat {
void registerCosseratTopologyBuilder(sofa::core::ObjectFactory *factory) {
    factory->registerObjects(
        sofa::core::ObjectRegistrationData("CosseratTopologyBuilder")
            .add<sofa::component::cosserat::engine::CosseratTopologyBuilder>());
}
}  // namespace Cosserat

namespace sofa::component::cosserat::engine {

CosseratTopologyBuilder::CosseratTopologyBuilder()
    : d_totalLength(initData(&d_totalLength, 15.0, "totalLength", "Total length of the beam")),
      d_nbSections(initData(&d_nbSections, 16, "nbSections", "Number of sections")),
      d_radius(initData(&d_radius, 0.1, "radius", "Radius of the cross section")),
      d_youngModulus(initData(&d_youngModulus, 1e6, "youngModulus", "Young's modulus (E)")),
      d_shearModulus(initData(&d_shearModulus, 1e5, "shearModulus", "Shear modulus (G)")),
      d_EA(initData(&d_EA, 0.0, "EA", "Axial stiffness (output)")),
      d_GA(initData(&d_GA, 0.0, "GA", "Shear stiffness (output)")),
      d_EIy(initData(&d_EIy, 0.0, "EIy", "Bending stiffness Y (output)")),
      d_EIz(initData(&d_EIz, 0.0, "EIz", "Bending stiffness Z (output)")),
      d_GJ(initData(&d_GJ, 0.0, "GJ", "Torsional stiffness (output)")),
      l_intrinsicState(initLink("intrinsicState", "Link to the CosseratIntrinsicState")),
      l_topology(initLink("topology", "Link to the EdgeSetTopologyContainer")) {}

void CosseratTopologyBuilder::init() {
    sofa::core::objectmodel::BaseObject::init();

    double length = d_totalLength.getValue();
    int sections = d_nbSections.getValue();
    double radius = d_radius.getValue();
    double E = d_youngModulus.getValue();
    double G = d_shearModulus.getValue();

    std::vector<sofa::type::Vec3d> nodes_vector;
    std::vector<SO3> frames_vector;
    std::vector<CosseratBeamGeometryEngine::Edge> edges_vector;
    std::vector<sofa::type::Vec3d> curvatures_vector(sections, sofa::type::Vec3d(0.0, 0.0, 0.0));

    // Call the geometry engine
    CosseratBeamGeometryEngine::buildConfiguration(length, sections, curvatures_vector,
                                                   nodes_vector, frames_vector, edges_vector);

    if (nodes_vector.size() != frames_vector.size() + 1) {
        msg_error() << "Geometry inconsistency: expected Size(X) = Size(R) + 1, got "
                    << nodes_vector.size() << " vs " << frames_vector.size();
        return;
    }

    // Calculate stiffness parameters
    double ea, ga, eiy, eiz, gj;
    CosseratBeamGeometryEngine::calculateStiffnessMatrix(E, G, radius, ea, ga, eiy, eiz, gj);

    d_EA.setValue(ea);
    d_GA.setValue(ga);
    d_EIy.setValue(eiy);
    d_EIz.setValue(eiz);
    d_GJ.setValue(gj);

    // Transfer positions and orientations to CosseratIntrinsicState
    if (l_intrinsicState.get()) {
        sofa::type::vector<sofa::type::Vec3d> positions(nodes_vector.begin(), nodes_vector.end());
        // d_orientations stores rotation vectors ω = log(R) ∈ ℝ³ (SOFA-accessible from Python)
        sofa::type::vector<sofa::type::Vec3d> orientations;
        orientations.reserve(frames_vector.size());
        for (const auto& R : frames_vector) {
            const auto omega = R.log();
            orientations.emplace_back(omega.x(), omega.y(), omega.z());
        }
        l_intrinsicState.get()->d_positions.setValue(positions);
        l_intrinsicState.get()->d_orientations.setValue(orientations);
    } else {
        msg_warning() << "CosseratIntrinsicState link not set — skipping DOF initialisation.";
    }

    if (l_topology.get()) {
        // Build edge connectivity as SOFA topology edges (unsigned int pairs)
        using TopoEdge = sofa::component::topology::container::dynamic::EdgeSetTopologyContainer::Edge;
        sofa::type::vector<TopoEdge> topo_edges;
        topo_edges.reserve(edges_vector.size());
        for (const auto& e : edges_vector) {
            topo_edges.push_back(TopoEdge{e.first, e.second});
        }
        l_topology.get()->d_edge.setValue(topo_edges);
    } else {
        msg_warning() << "EdgeSetTopologyContainer link not set — skipping edge initialisation.";
    }
}

}  // namespace sofa::component::cosserat::engine
