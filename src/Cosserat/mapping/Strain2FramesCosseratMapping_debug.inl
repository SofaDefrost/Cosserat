/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                                                                             *
 * This file is part of the Cosserat plugin.                                   *
 ******************************************************************************/
#pragma once

/**
 * @file Strain2FramesCosseratMapping_debug.inl
 * @brief Out-of-line implementations of the display*() debug helpers.
 *
 * These ~250 lines used to live at the bottom of Strain2FramesCosseratMapping.inl
 * but cluttered the file without adding to the mapping's correctness or runtime
 * behaviour.  They are only invoked when `d_debug.getValue() == true` (or via
 * `f_printLog`), so the cost in release builds is one branch + the compiled
 * std::cout instructions — kept here so the main .inl stays focused on the
 * mapping algebra.
 *
 * Included unconditionally at the end of Strain2FramesCosseratMapping.inl.
 * The implementations remain in a header because they're class-template
 * member functions and must be visible to every translation unit that
 * instantiates Strain2FramesCosseratMapping.
 */

#include <Cosserat/mapping/Strain2FramesCosseratMapping.h>

namespace Cosserat::mapping {

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displayStrainState(
			const sofa::type::vector<Coord1> &strainState,
			const std::string &context) const {
		std::cout << "\n=== STRAIN STATE DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Strain state size: " << strainState.size() << std::endl;

		for (size_t i = 0; i < strainState.size(); ++i) {
			std::cout << "  Strain[" << i << "]: ";

			if constexpr (std::is_same_v<Coord1, sofa::type::Vec3>) {
				std::cout << "[" << strainState[i][0] << ", " << strainState[i][1] << ", " << strainState[i][2] << "]";
			} else {
				std::cout << "[";
				for (int j = 0; j < strainState[i].size(); ++j) {
					std::cout << strainState[i][j];
					if (j < strainState[i].size() - 1)
						std::cout << ", ";
				}
				std::cout << "]";
			}
			std::cout << std::endl;
		}
		std::cout << "================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displayRigidState(
			const sofa::type::vector<sofa::Coord_t<In2>> &rigidState,
			const std::string &context) const {
		std::cout << "\n=== RIGID STATE DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Rigid state size: " << rigidState.size() << std::endl;

		for (size_t i = 0; i < rigidState.size(); ++i) {
			const auto &coord = rigidState[i];
			const auto &center = coord.getCenter();
			const auto &orientation = coord.getOrientation();

			std::cout << "  Rigid[" << i << "]:";
			std::cout << " pos=[" << center[0] << ", " << center[1] << ", " << center[2] << "]";
			std::cout << " quat=[" << orientation[0] << ", " << orientation[1] << ", " << orientation[2] << ", "
					  << orientation[3] << "]" << std::endl;
		}
		std::cout << "==============================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displayOutputFrames(
			const sofa::type::vector<OutCoord> &outputFrames,
			const std::string &context) const {
		std::cout << "\n=== OUTPUT FRAMES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Output frames size: " << outputFrames.size() << std::endl;

		for (size_t i = 0; i < outputFrames.size(); ++i) {
			const auto &frame = outputFrames[i];
			const auto &center = frame.getCenter();
			const auto &orientation = frame.getOrientation();

			std::cout << "  Frame[" << i << "]:";
			std::cout << " pos=[" << center[0] << ", " << center[1] << ", " << center[2] << "]";
			std::cout << " quat=[" << orientation[0] << ", " << orientation[1] << ", " << orientation[2] << ", "
					  << orientation[3] << "]" << std::endl;

			if (i > 0) {
				sofa::type::Vec3 diff = center - outputFrames[i - 1].getCenter();
				std::cout << "    Distance to prev: " << diff.norm() << std::endl;
			}
		}
		std::cout << "==================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displaySectionProperties(
			const std::string &context) const {
		std::cout << "\n=== SECTION PROPERTIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Section properties size: " << m_section_properties.size() << std::endl;

		for (size_t i = 0; i < m_section_properties.size(); ++i) {
			const auto &section = m_section_properties[i];
			const auto &strain = section.getStrainsVec();
			const auto &transform = section.getTransformation();

			std::cout << "  Section[" << i << "]:";
			std::cout << " length=" << section.getLength();
			std::cout << " strain=[" << strain << "]";
			std::cout << " indices=[" << section.getIndex0() << ", " << section.getIndex1() << "]" << std::endl;

			const auto &translation = transform.translation();
			const auto &rotation = transform.rotation();
			std::cout << "    Transform: trans=[" << translation[0] << ", " << translation[1] << ", " << translation[2]
					  << "]";
			std::cout << " rot_det=" << rotation.matrix().determinant() << std::endl;
		}
		std::cout << "=====================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displayFrameProperties(
			const std::string &context) const {
		std::cout << "\n=== FRAME PROPERTIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Frame properties size: " << m_frame_properties.size() << std::endl;

		for (size_t i = 0; i < m_frame_properties.size(); ++i) {
			const auto &frame = m_frame_properties[i];
			const auto &transform = frame.getTransformation();

			std::cout << "  Frame[" << i << "]:";
			std::cout << " length=" << frame.getLength();
			std::cout << " frames_sect_length_=" << frame.getLength();

			if (i < m_frame_to_section_indices.size()) {
				std::cout << " section_idx=" << m_frame_to_section_indices[i];
			}

			std::cout << " distance_to_nearest_beam_node=" << frame.getDistanceToSectionStart();

			const auto &translation = transform.translation();
			const auto &rotation = transform.rotation();
			std::cout << " trans=[" << translation[0] << ", " << translation[1] << ", " << translation[2] << "]";
			std::cout << " rot_det=" << rotation.matrix().determinant() << std::endl;
		}
		std::cout << "===================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displaySE3Transform(
			const SE3Types &transform,
			const std::string &name) const {
		std::cout << "\n=== SE3 TRANSFORM DEBUG: " << name << " ===\n";

		const auto &translation = transform.translation();
		const auto &rotation = transform.rotation();

		std::cout << "Translation: [" << translation[0] << ", " << translation[1] << ", " << translation[2] << "]\n";
		std::cout << "Rotation matrix:\n";
		const auto &R = rotation.matrix();
		for (int i = 0; i < 3; ++i) {
			std::cout << "  [" << R(i, 0) << ", " << R(i, 1) << ", " << R(i, 2) << "]\n";
		}
		std::cout << "Rotation determinant: " << R.determinant() << std::endl;

		Eigen::Quaternion<double> quat(R);
		std::cout << "Quaternion: [" << quat.w() << ", " << quat.x() << ", " << quat.y() << ", " << quat.z() << "]\n";

		std::cout << "Matrix form:\n";
		const auto &M = transform.matrix();
		for (int i = 0; i < 4; ++i) {
			std::cout << "  [" << M(i, 0) << ", " << M(i, 1) << ", " << M(i, 2) << ", " << M(i, 3) << "]\n";
		}
		std::cout << "==========================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displayMappingState(
			const std::string &context) const {
		std::cout << "\n=== MAPPING STATE DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Base index: " << d_baseIndex.getValue() << std::endl;
		std::cout << "Debug mode: " << (d_debug.getValue() ? "ON" : "OFF") << std::endl;

		// TODO: switch to m_frame_properties once frame topology is exposed there.
		const auto &curvFrames = d_curv_abs_frames.getValue();
		std::cout << "Curv abs frames size: " << curvFrames.size() << std::endl;
		if (!curvFrames.empty()) {
			std::cout << "  Values: [";
			for (size_t i = 0; i < curvFrames.size(); ++i) {
				std::cout << curvFrames[i];
				if (i < curvFrames.size() - 1)
					std::cout << ", ";
			}
			std::cout << "]\n";
		}

		std::cout << "Indices vectors size: " << m_frame_to_section_indices.size() << std::endl;
		if (!m_frame_to_section_indices.empty()) {
			std::cout << "  Values: [";
			for (size_t i = 0; i < m_frame_to_section_indices.size(); ++i) {
				std::cout << m_frame_to_section_indices[i];
				if (i < m_frame_to_section_indices.size() - 1)
					std::cout << ", ";
			}
			std::cout << "]\n";
		}

		std::cout << "Beam length vectors size: " << m_section_length_vectors.size() << std::endl;
		if (!m_section_length_vectors.empty()) {
			std::cout << "  Values: [";
			for (size_t i = 0; i < m_section_length_vectors.size(); ++i) {
				std::cout << m_section_length_vectors[i];
				if (i < m_section_length_vectors.size() - 1)
					std::cout << ", ";
			}
			std::cout << "]\n";
		}

		std::cout << "==============================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::displayVelocities(
			const sofa::type::vector<Deriv1> &strainVel,
			const sofa::type::vector<sofa::Deriv_t<In2>> &baseVel,
			const sofa::type::vector<OutDeriv> &outputVel,
			const std::string &context) const {
		std::cout << "\n=== VELOCITIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";

		std::cout << "Strain velocities (size: " << strainVel.size() << "):" << std::endl;
		for (size_t i = 0; i < strainVel.size(); ++i) {
			std::cout << "  StrainVel[" << i << "]: ";

			if constexpr (std::is_same_v<Deriv1, sofa::type::Vec3>) {
				std::cout << "[" << strainVel[i][0] << ", " << strainVel[i][1] << ", " << strainVel[i][2] << "]";
			} else {
				std::cout << "[";
				for (int j = 0; j < strainVel[i].size(); ++j) {
					std::cout << strainVel[i][j];
					if (j < strainVel[i].size() - 1)
						std::cout << ", ";
				}
				std::cout << "]";
			}
			std::cout << std::endl;
		}

		std::cout << "\nBase velocities (size: " << baseVel.size() << "):" << std::endl;
		for (size_t i = 0; i < baseVel.size(); ++i) {
			const auto &vel = baseVel[i];
			std::cout << "  BaseVel[" << i << "]: [" << vel[0] << ", " << vel[1] << ", " << vel[2] << ", " << vel[3]
					  << ", " << vel[4] << ", " << vel[5] << "]" << std::endl;
		}

		std::cout << "\nOutput velocities (size: " << outputVel.size() << "):" << std::endl;
		for (size_t i = 0; i < outputVel.size(); ++i) {
			const auto &vel = outputVel[i];
			std::cout << "  OutputVel[" << i << "]: [" << vel[0] << ", " << vel[1] << ", " << vel[2] << ", " << vel[3]
					  << ", " << vel[4] << ", " << vel[5] << "]" << std::endl;
		}

		std::cout << "==========================\n";
	}

} // namespace Cosserat::mapping
