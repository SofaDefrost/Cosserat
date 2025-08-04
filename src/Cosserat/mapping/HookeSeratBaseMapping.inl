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

			// Frame info initialized

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
			node.setStrain(strain[i]);
			m_section_properties.push_back(node);
		}
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
			}
		}
	}

	template<class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::computeTangExpImplementation(const double& curv_abs,
	const TangentVector & strain, const AdjointMatrix &adjoint_matrix, AdjointMatrix & tang_adjoint_matrix)
	{
		SReal theta = Vector3(strain(0), strain(1), strain(2)).norm();
		//old method
		//Matrix3 tilde_k = SO3Type::buildAntisymmetric(Vector3(strain_i(0), strain_i(1), strain_i(2)));// getTildeMatrix(Vec3(strain_i(0), strain_i(1), strain_i(2)));
		//Matrix3 tilde_q = SO3Type::buildAntisymmetric(Vector3(strain_i(3), strain_i(4), strain_i(5)));
		//buildAdjoint(tilde_k, tilde_q, ad_Xi);

		//@info: new method with Lie algebra
		//AdjointMatrix adjoint_matrix = node_info.getAdjoint();


		//@todo : compare the result of these two methods for computing
		//@todo : adjoint matrix;

		tang_adjoint_matrix = AdjointMatrix::Zero();
		AdjointMatrix Id6 = AdjointMatrix::Identity();

		if (theta <= std::numeric_limits<double>::epsilon()) {
			double scalar0 = std::pow(curv_abs, 2) / 2.0;
			tang_adjoint_matrix = curv_abs * Id6 + scalar0 * adjoint_matrix;
		} else {
			double scalar1 = (4.0 - 4.0 * cos(curv_abs * theta) -
							  curv_abs * theta * sin(curv_abs * theta)) /
							 (2.0 * theta * theta);
			double scalar2 = (4.0 * curv_abs * theta +
							  curv_abs * theta * cos(curv_abs * theta) -
							  5.0 * sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta);
			double scalar3 = (2.0 - 2.0 * cos(curv_abs * theta) -
							  curv_abs * theta * sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta * theta);
			double scalar4 = (2.0 * curv_abs * theta +
							  curv_abs * theta * cos(curv_abs * theta) -
							  3.0 * sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta * theta * theta);

			tang_adjoint_matrix = curv_abs * Id6 + scalar1 * adjoint_matrix + scalar2 * adjoint_matrix * adjoint_matrix +
				  scalar3 * adjoint_matrix * adjoint_matrix * adjoint_matrix +
				  scalar4 * adjoint_matrix * adjoint_matrix * adjoint_matrix * adjoint_matrix;
		}
	}


	template<class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::updateTangExpSE3() {

		auto node_count = m_section_properties.size();

		//update node's tang SE3 matrix
		AdjointMatrix tang_matrix = AdjointMatrix::Zero();
		m_section_properties[0].setTanAdjointMatrix(tang_matrix);

		for (auto i = 1; i< node_count; i++ ) {
			auto node_info = m_section_properties[i];
			computeTangExpImplementation(
				node_info.getLength(),
				node_info.getStrainsVec(),
				node_info.getAdjoint(),
				tang_matrix);
			node_info.setTanAdjointMatrix(tang_matrix);
			std::cout << "Node[" << i << "] tang adjoint matrix: \n" << tang_matrix << std::endl;
		}

		//update frames's tang SE3 matrix
		auto frame_count = m_frameProperties.size();
		for (auto i = 0; i<frame_count; i++) {
			auto frame_info = m_frameProperties[i];
			auto related_section_index = m_frameProperties[i].get_related_beam_index_();
			auto frame_strain = m_section_properties[related_section_index].getStrainsVec();
			computeTangExpImplementation(
				frame_info.getDistanceToNearestBeamNode(),
				frame_strain,
				frame_info.getAdjoint(),
				tang_matrix);
			frame_info.setTanAdjointMatrix(tang_matrix);
			std::cout << "Frame[" << i << "] tang adjoint matrix: \n"
			<< tang_matrix << std::endl;
		}


	}
} // namespace Cosserat::mapping
