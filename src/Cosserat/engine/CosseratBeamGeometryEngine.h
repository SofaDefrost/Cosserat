#pragma once

#include <liegroups/SO3.h>
#include <sofa/type/Vec.h>

#include <utility>
#include <vector>

namespace sofa::component::cosserat::engine {

// 1. Structure for geometry parameters (equivalent to GeometryParams in Python)
struct CosseratBeamGeometryParameters {
    double radius = 0.1;
    int nbSections = 16;
    int nbFrames = 15;
    double totalLength = 15.0;
};

// 2. CosseratBeamGeometryEngine Class
class CosseratBeamGeometryEngine {
   public:
    using Vec3d = sofa::type::Vec3d;
    using SO3 = sofa::component::cosserat::liegroups::SO3<double>;
    using Edge = std::pair<unsigned int, unsigned int>;

    /**
     * @brief Builds a straight configuration of the beam along the X axis.
     *
     * @param totalLength Total length of the beam.
     * @param nbSections Number of sections (segments).
     * @param nodes Output vector of node positions.
     * @param frames Output vector of node orientations (SO3).
     * @param edges Output vector of edges (connectivity).
     */
    static void buildStraightConfiguration(double totalLength, int nbSections,
                                           std::vector<Vec3d> &nodes, std::vector<SO3> &frames,
                                           std::vector<Edge> &edges);

    /**
     * @brief Builds a curved configuration of the beam using the Lie exponential.
     *
     * @param totalLength Total length of the beam.
     * @param nbSections Number of sections (segments).
     * @param curvatures Curvature vector (Omega = [torsion_x, bending_y, bending_z]) for each segment.
     * @param nodes Output vector of node positions.
     * @param frames Output vector of node orientations (SO3).
     * @param edges Output vector of edges (connectivity).
     */
    static void buildCurvedConfiguration(double totalLength, int nbSections,
                                         const std::vector<Vec3d> &curvatures,
                                         std::vector<Vec3d> &nodes, std::vector<SO3> &frames,
                                         std::vector<Edge> &edges);

    /**
     * @brief Builds the intrinsic configuration with N+1 nodes and N rotations.
     *
     * @param totalLength Total length of the beam.
     * @param nbSections Number of sections (segments).
     * @param curvatures Curvature vector (Omega = [torsion_x, bending_y, bending_z]) for each segment.
     * @param nodes Output vector of node positions (size N+1).
     * @param frames Output vector of segment orientations (size N).
     * @param edges Output vector of edges (connectivity).
     */
    static void buildConfiguration(double totalLength, int nbSections,
                                   const std::vector<Vec3d> &curvatures, std::vector<Vec3d> &nodes,
                                   std::vector<SO3> &frames, std::vector<Edge> &edges);

    /**
     * @brief Calculates the stiffness matrix constants for a circular cross-section.
     *
     * @param E Young's modulus.
     * @param G Shear modulus.
     * @param r Radius of the section.
     * @param EA Axial stiffness.
     * @param GA Shear stiffness.
     * @param EIy Bending stiffness (Y axis).
     * @param EIz Bending stiffness (Z axis).
     * @param GJ Torsional stiffness.
     */
    static void calculateStiffnessMatrix(double E, double G, double r, double &EA, double &GA,
                                         double &EIy, double &EIz, double &GJ);

    /**
     * @brief Calculates the stiffness matrix constants for a rectangular cross-section.
     *
     * @param E Young's modulus.
     * @param G Shear modulus.
     * @param w Width of the section.
     * @param h Height of the section.
     * @param EA Axial stiffness.
     * @param GA Shear stiffness.
     * @param EIy Bending stiffness (Y axis).
     * @param EIz Bending stiffness (Z axis).
     * @param GJ Torsional stiffness.
     */
    static void calculateStiffnessMatrix(double E, double G, double w, double h, double &EA,
                                         double &GA, double &EIy, double &EIz, double &GJ);
};

}  // namespace sofa::component::cosserat::engine
