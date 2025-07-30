/*
 * HookeSeratBaseMapping.inl
 * Implementation details for the HookeSeratBaseMapping class.
 * This file contains implementations for functions inline in the HookeSeratBaseMapping class.
 */
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <iostream>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>

namespace Cosserat::mapping {

	using namespace sofa::component::cosserat::liegroups;
	using sofa::helper::getReadAccessor;
	using sofa::type::Quat;
	using sofa::type::Vec3;
	using sofa::type::Vec6;
	using sofa::type::vector;

	template<class TIn1, class TIn2, class TOut>
	HookeSeratBaseMapping<TIn1, TIn2, TOut>::HookeSeratBaseMapping() :
		d_curv_abs_section(initData(&d_curv_abs_section, "curv_abs_input",
									"Curvilinear abscissa of the input sections along the rod")),
		d_curv_abs_frames(initData(&d_curv_abs_frames, "curv_abs_output",
								   "Curvilinear abscissa of the output frames along the rod")),
		d_debug(initData(&d_debug, false, "debug", "Enable debug output")), m_strain_state(nullptr),
		m_rigid_base(nullptr), m_frames(nullptr) {
		msg_info("HookeSeratBaseMapping") << "HookeSeratBaseMapping constructor called  !!!";
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::init() {
		// Initialize pointers to nullptr
		msg_info("HookeSeratBaseMapping") << "Initializing HookeSeratBaseMapping...";

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

		// Get the initial configuration g(X):frames and initialize FrameInfo objects
		if (m_frames) {
			auto xfromData = m_frames->read(sofa::core::vec_id::read_access::position);
			const auto &xfrom = xfromData->getValue();

			// Initialize frame properties using the initial frame states
			const auto frame_count = xfrom.size();

			m_frameProperties.clear();
			m_frameProperties.reserve(frame_count);

			for (size_t i = 0; i < frame_count; ++i) {
				// Convert SOFA coordinates to SE3 transformations
				const auto &frame_i = xfrom[i];
				Vector3 translation(frame_i.getCenter()[0], frame_i.getCenter()[1], frame_i.getCenter()[2]);

				// Convert quaternion to rotation matrix
				const auto &quat = frame_i.getOrientation();
				SO3Type rotation;
				// Convert SOFA quaternion to our SO3 representation
				// SOFA quaternions use [x, y, z, w] order, Eigen uses [w, x, y, z]
				Eigen::Quaternion<double> eigenQuat(quat[3], quat[0], quat[1], quat[2]);
				rotation = SO3Type(eigenQuat.toRotationMatrix());

				SE3Type gXi(rotation, translation);

				std::cout << " frame i : "<< i << " gXi : "<< translation.transpose() << "  " << rotation<< std::endl;

				// Create FrameInfo with initial transformation
				// Length and kappa will be set later in initializeFrameProperties
				FrameInfo frameInfo;
				frameInfo.setTransformation(gXi);
				m_frameProperties.push_back(frameInfo);
			}
		}

		// Initialize geometry information
		updateGeometryInfo();

		// Initialize section and frame properties based on geometry
		initializeSectionProperties();

		// No need anymore, this is done in updateGeometryInfo
		initializeFrameProperties();

		// Validation
		if (!validateSectionProperties()) {
			msg_error() << "Invalid section properties detected";
			return;
		}

		// Check continuity (with warning only)
		// if (!checkContinuity()) {
		// 	msg_warning() << "Rod sections are not continuous";
		// }

		// Pre-compute adjoint matrices for performance
		for (const auto &section: m_section_properties) {
			section.getAdjoint(); // Force computation and caching
		}

		// Call parent initialization
		Inherit::init();
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::updateGeometryInfo() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto &curv_abs_frames = d_curv_abs_frames.getValue();

		if (curv_abs_frames.empty()) {
			msg_warning("HookeSeratBaseMapping") << "Empty frames data";
			return;
		}

		const auto frame_count = curv_abs_frames.size();

		// Pré-allocation efficace
		reserveContainers(frame_count);

		size_t current_section_index = 1;
		constexpr double TOLERANCE = 1e-3;

		for (size_t i = 0; i < frame_count; ++i) {
			const auto frame_pos = curv_abs_frames[i];

			// Find the index of the section of on which each frame belong on
			auto result = findSectionIndex(frame_pos, curv_abs_section, current_section_index, TOLERANCE);
			// update the beam's frames strcut
			updateFrameData(i, result.index_for_frame, frame_pos, curv_abs_section);

			// Mettre à jour pour la prochaine itération
			current_section_index = result.index_for_next;
		}
		logCompletionInfo();
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::initializeSectionProperties() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto node_count = curv_abs_section.size();
		const auto &strain = m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

		// Initialize section properties based on geometry
		m_section_properties.clear();
		m_section_properties.reserve(node_count);

		// Save node 0 info
		TangentVector init_strain = TangentVector::Zero(); // Initial curvature, will be updated
		SectionInfo node_0(0., init_strain, 0, SE3Type::computeIdentity()); // First node has zero length
		// The first node is often fix or attached to an object
		m_section_properties.push_back(node_0);

		// Initial strain vector, can be modified later
		TangentVector strain_0 = TangentVector::Zero();

		// compute the length of each beam segment.
		//Initialize a 6D strain vector (angular and linear components)
		std::adjacent_difference(curv_abs_section.begin() + 1, curv_abs_section.end(),
								 std::back_inserter(m_beam_length_vectors));

		for (size_t i = 0; i < node_count - 1; ++i) {
			double length = m_beam_length_vectors[i];
			// Fill with

			SectionInfo node(length, strain_0, i, SE3Type::computeIdentity());
			std::cout << "===> strain[i] : " << strain[i] << std::endl;
			node.setStrain(strain[i]);
			m_section_properties.push_back(node);
		}
		std::cout << "HookeSeratBaseMapping" << "  m_beam_length_vectors: " << m_beam_length_vectors.size()
				  << " elements" << std::endl;
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::initializeFrameProperties() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto &curv_abs_frames = d_curv_abs_frames.getValue();

		// Initialize frame properties based on the number of frames
		// The initiation is done in previous updateGeometryInfo

		// Create frame properties for each frame
		for (size_t i = 0; i < curv_abs_frames.size()-1; ++i) {
			if (i < m_indices_vectors.size()) {
				// Calculate frame section length based on position relative to beam nodes
				// This should be the distance from the frame to the closest beam node toward the base
				double frame_length = curv_abs_frames[i + 1] - curv_abs_frames[i];

				// Ensure frameLength is positive (safety check)
				if (frame_length <= 0) {
					msg_warning() << "Frame " << i << " has non-positive length " << frame_length
								  << ". Frame pos: " << curv_abs_frames[i]
								  << ", Section pos: " << curv_abs_section[m_indices_vectors[i] - 1]
								  << ". Using curv_abs_frames[i] instead.";
					frame_length = curv_abs_frames[i];
				}

				m_frameProperties[i+1].setLength(frame_length);
				std::cout<< "Frame : "<< i+1 << "  length : " << frame_length << std::endl;
			}
		}
	}

} // namespace Cosserat::mapping
