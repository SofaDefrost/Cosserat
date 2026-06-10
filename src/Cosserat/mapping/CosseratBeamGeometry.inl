/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                                                                             *
 * This file is part of the Cosserat plugin.                                   *
 ******************************************************************************/
#pragma once

#include <Cosserat/mapping/CosseratBeamGeometry.h>
#include <iterator>
#include <numeric>

namespace Cosserat::mapping {

	// ──────────────────────────────────────────────────────────────────────────
	// Geometry building blocks
	// ──────────────────────────────────────────────────────────────────────────

	template<class TStrain>
	void CosseratBeamGeometry<TStrain>::updateGeometryInfo() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto &curv_abs_frames  = d_curv_abs_frames.getValue();

		if (curv_abs_frames.empty()) {
			msg_warning("CosseratBeamGeometry") << "Empty frames data";
			return;
		}

		const auto frame_count = curv_abs_frames.size();
		reserveContainers(frame_count);

		size_t current_section_index = 1;
		constexpr double TOLERANCE = 1e-3;

		for (size_t i = 0; i < frame_count; ++i) {
			const auto frame_pos = curv_abs_frames[i];
			auto result = findSectionIndex(frame_pos, curv_abs_section, current_section_index, TOLERANCE);
			updateFrameData(i, result.index_for_frame, frame_pos, curv_abs_section);
			current_section_index = result.index_for_next;
		}
		logCompletionInfo();
	}

	template<class TStrain>
	void CosseratBeamGeometry<TStrain>::initializeSectionProperties(const StrainVecCoord &strain) {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto node_count = curv_abs_section.size();

		m_section_properties.clear();
		m_section_properties.reserve(node_count);

		// First node (anchor): zero length, identity transform
		TangentVector init_strain = TangentVector::Zero();
		SectionInfo node_0(0., init_strain, 0, SE3Type::computeIdentity());
		m_section_properties.push_back(node_0);

		TangentVector strain_0 = TangentVector::Zero();

		// Compute per-section lengths via adjacent_difference on curvilinear abscissas
		std::adjacent_difference(curv_abs_section.begin() + 1, curv_abs_section.end(),
								 std::back_inserter(m_section_length_vectors));

		for (size_t i = 0; i < node_count - 1; ++i) {
			double length = m_section_length_vectors[i];
			SectionInfo node(length, strain_0, i, SE3Type::computeIdentity());
			node.setIndices(i);
			node.setStrain(strain[i]);
			m_section_properties.push_back(node);
		}
	}

	template<class TStrain>
	void CosseratBeamGeometry<TStrain>::initializeFrameProperties() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto &curv_abs_frames  = d_curv_abs_frames.getValue();

		for (size_t i = 0; i + 1 < curv_abs_frames.size(); ++i) {
			if (i < m_frame_to_section_indices.size()) {
				double frame_length = curv_abs_frames[i + 1] - curv_abs_frames[i];

				if (frame_length <= 0) {
					msg_warning("CosseratBeamGeometry")
						<< "Frame " << i << " has non-positive length " << frame_length
						<< ". Frame pos: " << curv_abs_frames[i]
						<< ", Section pos: " << curv_abs_section[m_frame_to_section_indices[i] - 1]
						<< ". Using curv_abs_frames[i] instead.";
					frame_length = curv_abs_frames[i];
				}
				m_frame_properties[i + 1].setLength(frame_length);
			}
		}
	}

	template<class TStrain>
	void CosseratBeamGeometry<TStrain>::updateTangExpSE3() {
		auto node_count = m_section_properties.size();

		AdjointMatrix tang_matrix = AdjointMatrix::Zero();
		m_section_properties[0].setTanAdjointMatrix(tang_matrix);

		for (size_t i = 1; i < node_count; ++i) {
			auto &node_info = m_section_properties[i];
			computeTangExpImplementation(node_info.getLength(), node_info.getStrainsVec(), tang_matrix);
			node_info.setTanAdjointMatrix(tang_matrix);
			if (d_debug.getValue()) {
				msg_info("CosseratBeamGeometry") << "Node[" << i << "] tang adjoint matrix:\n" << tang_matrix;
			}
		}

		auto frame_count = m_frame_properties.size();
		for (size_t i = 0; i < frame_count; ++i) {
			auto &frame_info = m_frame_properties[i];
			auto related_section_index = frame_info.getRelatedSectionIndex();
			auto frame_strain = m_section_properties[related_section_index].getStrainsVec();
			computeTangExpImplementation(frame_info.getDistanceToSectionStart(), frame_strain, tang_matrix);
			frame_info.setTanAdjointMatrix(tang_matrix);
			if (d_debug.getValue()) {
				msg_info("CosseratBeamGeometry") << "Frame[" << i << "] tang adjoint matrix:\n" << tang_matrix;
			}
		}
	}

	// ──────────────────────────────────────────────────────────────────────────
	// Static helpers
	// ──────────────────────────────────────────────────────────────────────────

	template<class TStrain>
	Matrix3 CosseratBeamGeometry<TStrain>::getTildeMatrix(const Vector3 &u) {
		Matrix3 tild = Matrix3::Zero();
		tild(0, 1) = -u[2];
		tild(0, 2) =  u[1];
		tild(1, 2) = -u[0];
		tild(1, 0) = -tild(0, 1);
		tild(2, 0) = -tild(0, 2);
		tild(2, 1) = -tild(1, 2);
		return tild;
	}

	template<class TStrain>
	void CosseratBeamGeometry<TStrain>::buildAdjoint(const Matrix3 &A, const Matrix3 &B, AdjointMatrix &Adjoint) {
		Adjoint = AdjointMatrix::Zero();
		for (unsigned int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				Adjoint(i, j)         = A(i, j);
				Adjoint(i + 3, j + 3) = A(i, j);
				Adjoint(i + 3, j)     = B(i, j);
			}
		}
	}

	template<class TStrain>
	void CosseratBeamGeometry<TStrain>::computeTangExpImplementation(const double &curv_abs,
																	 const TangentVector &strain,
																	 AdjointMatrix &tang_adjoint_matrix) {
		SReal theta = Vector3(strain(0), strain(1), strain(2)).norm();

		Matrix3 tilde_k = getTildeMatrix(Vector3(strain(0), strain(1), strain(2)));
		Matrix3 tilde_q = getTildeMatrix(Vector3(strain(3), strain(4), strain(5)));
		AdjointMatrix ad_Xi = AdjointMatrix::Zero();
		buildAdjoint(tilde_k, tilde_q, ad_Xi);

		tang_adjoint_matrix = AdjointMatrix::Zero();
		AdjointMatrix Id6   = AdjointMatrix::Identity();

		if (theta <= std::numeric_limits<double>::epsilon()) {
			double scalar0 = std::pow(curv_abs, 2) / 2.0;
			tang_adjoint_matrix = curv_abs * Id6 + scalar0 * ad_Xi;
		} else {
			double scalar1 = (4.0 - 4.0 * std::cos(curv_abs * theta) -
							  curv_abs * theta * std::sin(curv_abs * theta)) /
							 (2.0 * theta * theta);
			double scalar2 = (4.0 * curv_abs * theta +
							  curv_abs * theta * std::cos(curv_abs * theta) -
							  5.0 * std::sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta);
			double scalar3 = (2.0 - 2.0 * std::cos(curv_abs * theta) -
							  curv_abs * theta * std::sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta * theta);
			double scalar4 = (2.0 * curv_abs * theta +
							  curv_abs * theta * std::cos(curv_abs * theta) -
							  3.0 * std::sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta * theta * theta);

			tang_adjoint_matrix = curv_abs * Id6 + scalar1 * ad_Xi + scalar2 * ad_Xi * ad_Xi +
								  scalar3 * ad_Xi * ad_Xi * ad_Xi +
								  scalar4 * ad_Xi * ad_Xi * ad_Xi * ad_Xi;
		}
	}

} // namespace Cosserat::mapping
