/*
 * HookeSeratBaseMapping.inl
 * Implementation details for the HookeSeratBaseMapping class.
 * This file contains implementations for functions inline in the HookeSeratBaseMapping class.
 */
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>

namespace Cosserat::mapping {

using namespace sofa::component::cosserat::liegroups;
using sofa::helper::getReadAccessor;
using sofa::type::Vec6;
using sofa::type::Vec3;
using sofa::type::Quat;
using sofa::type::vector;

template <class TIn1, class TIn2, class TOut>
HookeSeratBaseMapping<TIn1, TIn2, TOut>::HookeSeratBaseMapping()
    : d_curv_abs_section(initData(&d_curv_abs_section, "curv_abs_input",
                                  "Curvilinear abscissa of the input sections along the rod")),
      d_curv_abs_frames(initData(&d_curv_abs_frames, "curv_abs_output",
                                 "Curvilinear abscissa of the output frames along the rod")),
      d_debug(initData(&d_debug, false, "debug", "Enable debug output")),m_strain_state(nullptr),m_rigid_base(nullptr),m_frames(nullptr)
{
    // Initialize empty section and frame properties
    m_sectionProperties.clear();
    m_frameProperties.clear();
}

template <class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::init() {
    // Initialize pointers to nullptr
    m_strain_state = nullptr;
    m_rigid_base = nullptr;
    m_frames = nullptr;

    // Check if all required models are present
    if (this->fromModels1.empty()) {
        msg_error() << "Input1 (strain state) not found";
        return;
    }

    if (this->fromModels2.empty()) {
        msg_error() << "Input2 (rigid base) not found";
        return;
    }

    if (this->toModels.empty()) {
        msg_error() << "Output (frames) missing";
        return;
    }

	// Assign mechanical states
    m_strain_state = this->fromModels1[0];
    m_rigid_base = this->fromModels2[0];
    m_frames = this->toModels[0];

    // Get initial frame state and initialize FrameInfo objects
    if (m_frames) {
        auto xfromData = m_frames->read(sofa::core::vec_id::read_access::position);
        const auto& xfrom = xfromData->getValue();
        
        // Initialize frame properties using the initial frame states
        m_frameProperties.clear();
        m_frameProperties.reserve(xfrom.size());
        
        for (size_t i = 0; i < xfrom.size(); ++i) {
            // Convert SOFA coordinates to SE3 transformations
            const auto& coord = xfrom[i];
            Vector3 translation(coord.getCenter()[0], coord.getCenter()[1], coord.getCenter()[2]);
            
            // Convert quaternion to rotation matrix
            const auto& quat = coord.getOrientation();
            SO3Type rotation;
            // Convert SOFA quaternion to our SO3 representation
            // SOFA quaternions use [x, y, z, w] order, Eigen uses [w, x, y, z]
            Eigen::Quaternion<double> eigenQuat(quat[3], quat[0], quat[1], quat[2]);
            rotation = SO3Type(eigenQuat.toRotationMatrix());
            
            SE3Type transform(rotation, translation);
            
            // Create FrameInfo with initial transformation
            // Length and kappa will be set later in initializeFrameProperties
            FrameInfo frameInfo;
            frameInfo.setTransformation(transform);
            m_frameProperties.push_back(frameInfo);
        }
    }

    // Initialize geometry information
    updateGeometryInfo();
    
    // Initialize section and frame properties based on geometry
    initializeSectionProperties();
    initializeFrameProperties();

    // Validation
    if (!validateSectionProperties()) {
        msg_error() << "Invalid section properties detected";
        return;
    }

    // Check continuity (with warning only)
    if (!checkContinuity()) {
        msg_warning() << "Rod sections are not continuous";
    }

    // Pre-compute adjoint matrices for performance
    for (const auto& section : m_sectionProperties) {
        section.getAdjoint(); // Force computation and caching
    }

    // Call parent initialization
    Inherit::init();
}

template <class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::updateGeometryInfo() {
    // Similaire à BaseCosseratMapping::update_geometry_info() mais adapté pour les nouvelles structures
	const auto& curv_abs_section = d_curv_abs_section.getValue();
    const auto& curv_abs_frames = d_curv_abs_frames.getValue();

    msg_info() << "curv_abs_section: " << curv_abs_section.size()
               << "\ncurv_abs_frames: " << curv_abs_frames.size();

    // Clear existing geometry vectors
    m_indicesVectors.clear();
    m_indicesVectorsDraw.clear();
    m_beamLengthVectors.clear();

    const auto sz = curv_abs_frames.size();
    auto sectionIndex = 1;
    
    // Process frame indices similar to BaseCosseratMapping
    for (size_t i = 0; i < sz; ++i) {
        if (curv_abs_section[sectionIndex] > curv_abs_frames[i]) {
            m_indicesVectors.emplace_back(sectionIndex);
            m_indicesVectorsDraw.emplace_back(sectionIndex);
        }
        else if (std::abs(curv_abs_section[sectionIndex] - curv_abs_frames[i]) < 1e-6) {
            m_indicesVectors.emplace_back(sectionIndex);
            sectionIndex++;
            m_indicesVectorsDraw.emplace_back(sectionIndex);
        }
        else {
            sectionIndex++;
            m_indicesVectors.emplace_back(sectionIndex);
            m_indicesVectorsDraw.emplace_back(sectionIndex);
        }
    }

    // Calculate beam lengths
    for (size_t j = 0; j < sz - 1; ++j) {
        m_beamLengthVectors.emplace_back(curv_abs_section[j + 1] - curv_abs_section[j]);
    }

    msg_info() << "m_indicesVectors: " << m_indicesVectors.size() << " elements";
}

template <class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::initializeSectionProperties() {
    const auto& curv_abs_section = d_curv_abs_section.getValue();
    
    // Initialize section properties based on geometry
    m_sectionProperties.clear();
    m_sectionProperties.reserve(m_beamLengthVectors.size());
    
    for (size_t i = 0; i < m_beamLengthVectors.size(); ++i) {
        double length = m_beamLengthVectors[i];
        Vector3 kappa = Vector3::Zero(); // Initial curvature, will be updated
        
        SectionInfo section(length, kappa, i, i + 1);
        m_sectionProperties.push_back(section);
    }
}

template <class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::initializeFrameProperties() {
    const auto& curv_abs_section = d_curv_abs_section.getValue();
    const auto& curv_abs_frames = d_curv_abs_frames.getValue();
    
    // Update frame properties with correct lengths and indices
    for (size_t i = 0; i < m_frameProperties.size(); ++i) {
        if (i < m_indicesVectors.size()) {
            // Calculate frame section length based on position relative to beam nodes
            double frameLength = curv_abs_frames[i] - curv_abs_section[m_indicesVectors[i] - 1];
            m_frameProperties[i].setLength(frameLength);
            
            // Set initial kappa (will be updated during simulation)
            Vector6 kappa = Vector6::Zero();
            m_frameProperties[i].setKappa(kappa);
        }
    }
}

} // namespace Cosserat::mapping

