#include "CosseratBeamGeometryEngine.h"

#include <algorithm>
#include <cmath>

namespace sofa::component::cosserat::engine {

void CosseratBeamGeometryEngine::buildStraightConfiguration(double totalLength, int nbSections,
                                                            std::vector<Vec3d> &nodes,
                                                            std::vector<SO3> &frames,
                                                            std::vector<Edge> &edges) {
    nodes.clear();
    frames.clear();
    edges.clear();

    if (nbSections <= 0 || totalLength <= 0.0) return;

    double sectionLength = totalLength / nbSections;
    int nbNodes = nbSections + 1;

    nodes.reserve(nbNodes);
    frames.reserve(nbSections);
    edges.reserve(nbSections);

    SO3 identity = SO3::identity();

    for (int i = 0; i < nbNodes; ++i) {
        nodes.emplace_back(Vec3d(i * sectionLength, 0.0, 0.0));
        if (i < nbSections) {
            frames.emplace_back(identity);
            edges.emplace_back(i, i + 1);
        }
    }
}

void CosseratBeamGeometryEngine::buildCurvedConfiguration(double totalLength, int nbSections,
                                                          const std::vector<Vec3d> &curvatures,
                                                          std::vector<Vec3d> &nodes,
                                                          std::vector<SO3> &frames,
                                                          std::vector<Edge> &edges) {
    nodes.clear();
    frames.clear();
    edges.clear();

    if (nbSections <= 0 || totalLength <= 0.0) return;

    double sectionLength = totalLength / nbSections;
    int nbNodes = nbSections + 1;

    nodes.reserve(nbNodes);
    frames.reserve(nbNodes);
    edges.reserve(nbSections);

    // Initialize R0 and x0
    nodes.emplace_back(Vec3d(0.0, 0.0, 0.0));
    frames.emplace_back(SO3::identity());

    for (int i = 0; i < nbSections; ++i) {
        // Retrieve the curvature for the current segment (or 0 by default)
        Vec3d omega(0.0, 0.0, 0.0);
        if (i < static_cast<int>(curvatures.size())) {
            omega = curvatures[i];
        }

        SO3::TangentVector omega_tangent(omega.x(), omega.y(), omega.z());

        SO3 R_i = frames[i];

        // R_{i+1} = R_i * exp(Omega_i)
        SO3 R_next = R_i * SO3::exp(omega_tangent);

        // x_{i+1} = x_i + R_i * [L_seg, 0, 0]^T
        SO3::Vector step_local(sectionLength, 0.0, 0.0);
        SO3::Vector step_global = R_i.act(step_local);

        Vec3d next_pos = nodes[i] + Vec3d(step_global.x(), step_global.y(), step_global.z());

        nodes.emplace_back(next_pos);
        frames.emplace_back(R_next);
        edges.emplace_back(i, i + 1);
    }
}

void CosseratBeamGeometryEngine::buildConfiguration(double totalLength, int nbSections,
                                                    const std::vector<Vec3d> &curvatures,
                                                    std::vector<Vec3d> &nodes,
                                                    std::vector<SO3> &frames,
                                                    std::vector<Edge> &edges) {
    nodes.clear();
    frames.clear();
    edges.clear();

    if (nbSections <= 0 || totalLength <= 0.0) return;

    double sectionLength = totalLength / nbSections;
    int nbNodes = nbSections + 1;

    nodes.reserve(nbNodes);
    frames.reserve(nbSections);
    edges.reserve(nbSections);

    // Initialize x0 and the temporary reference frame
    nodes.emplace_back(Vec3d(0.0, 0.0, 0.0));
    SO3 current_frame = SO3::identity();

    for (int i = 0; i < nbSections; ++i) {
        // Retrieve the curvature for the current segment (or 0 by default)
        Vec3d omega(0.0, 0.0, 0.0);
        if (i < static_cast<int>(curvatures.size())) {
            omega = curvatures[i];
        }

        SO3::TangentVector omega_tangent(omega.x(), omega.y(), omega.z());

        // R_i = R_{i-1} * exp(Omega_i)
        current_frame = current_frame * SO3::exp(omega_tangent);

        // Store the frame of the segment
        frames.emplace_back(current_frame);

        // x_{i+1} = x_i + R_i * [L_seg, 0, 0]^T
        SO3::Vector step_local(sectionLength, 0.0, 0.0);
        SO3::Vector step_global = current_frame.act(step_local);

        Vec3d next_pos = nodes[i] + Vec3d(step_global.x(), step_global.y(), step_global.z());

        nodes.emplace_back(next_pos);
        edges.emplace_back(i, i + 1);
    }
}

void CosseratBeamGeometryEngine::calculateStiffnessMatrix(double E, double G, double r, double &EA,
                                                          double &GA, double &EIy, double &EIz,
                                                          double &GJ) {
    double area = M_PI * r * r;
    double Iy = M_PI * std::pow(r, 4) / 4.0;
    double Iz = M_PI * std::pow(r, 4) / 4.0;
    double J = M_PI * std::pow(r, 4) / 2.0;

    EA = E * area;
    GA = G * area;
    EIy = E * Iy;
    EIz = E * Iz;
    GJ = G * J;
}

void CosseratBeamGeometryEngine::calculateStiffnessMatrix(double E, double G, double w, double h,
                                                          double &EA, double &GA, double &EIy,
                                                          double &EIz, double &GJ) {
    double area = w * h;
    double Iy = w * std::pow(h, 3) / 12.0;
    double Iz = h * std::pow(w, 3) / 12.0;

    // Approximation of the torsional constant for a rectangular cross-section
    double a = std::max(w, h);
    double b = std::min(w, h);
    double J =
        a * std::pow(b, 3) * (1.0 / 3.0 - 0.21 * (b / a) * (1.0 - std::pow(b / a, 4) / 12.0));

    EA = E * area;
    GA = G * area;
    EIy = E * Iy;
    EIz = E * Iz;
    GJ = G * J;
}

}  // namespace sofa::component::cosserat::engine
