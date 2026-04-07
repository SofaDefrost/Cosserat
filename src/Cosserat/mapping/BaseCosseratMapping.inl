/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.        *
 *******************************************************************************
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/BaseCosseratMapping.h>

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/Quat.h>

#include <string>

// To go further =>
// https://www.mathworks.com/matlabcentral/fileexchange/83038-sorosim/

namespace Cosserat::mapping {

	using sofa::helper::getReadAccessor;
	using sofa::type::Quat;
	using sofa::type::Vec3;
	using sofa::type::Vec6;

	template<class TIn1, class TIn2, class TOut>
	BaseCosseratMapping<TIn1, TIn2, TOut>::BaseCosseratMapping() :
		d_curv_abs_section(initData(&d_curv_abs_section, "curv_abs_input", " need to be com....")),
		d_curv_abs_frames(initData(&d_curv_abs_frames, "curv_abs_output", " need to be com....")),
		d_debug(initData(&d_debug, false, "debug", "printf for the debug")), m_index_input(0) {}


	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::init() {
		m_strain_state = nullptr;
		m_rigid_base = nullptr;
		m_global_frames = nullptr;

		if (fromModels1.empty()) {
			msg_error() << "input1 not found";
			return;
		}

		if (fromModels2.empty()) {
			msg_error() << "input2 not found";
			return;
		}

		if (toModels.empty()) {
			msg_error() << "output missing";
			return;
		}

		m_strain_state = fromModels1[0];
		m_rigid_base = fromModels2[0];
		m_global_frames = toModels[0];

		// Get initial frame state
		auto xfromData = m_global_frames->read(sofa::core::vec_id::read_access::position);
		const vector<OutCoord> xfrom = xfromData->getValue();

		m_vec_transform.clear();
		for (unsigned int i = 0; i < xfrom.size(); i++)
			m_vec_transform.push_back(xfrom[i]);

		// update_geometry_info();
		doBaseCosseratInit();
		// Inherit1::init();
	}

	//___________________________________________________________________________
	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::update_geometry_info() {
		// For each frame in the global frame, find the segment of the beam to which
		// it is attached. Here we only use the information from the curvilinear
		// abscissa of each frame.
		auto curv_abs_section = getReadAccessor(d_curv_abs_section);
		auto curv_abs_frames = getReadAccessor(d_curv_abs_frames);
		if (curv_abs_section.empty() || curv_abs_frames.empty()) {
			msg_warning() << "Empty curvilinear abscissa data";
			return;
		}
		if (curv_abs_section.size() < 2) {
			msg_error() << "Need at least 2 sections for beam geometry";
			return;
		}

		msg_info() << " curv_abs_section: " << curv_abs_section.size()
				   << "\ncurv_abs_frames: " << curv_abs_frames.size();

		std::cout << "==> Curv abs frames: " << curv_abs_frames << std::endl;
		std::cout << "==> Strain state: " << curv_abs_section << std::endl;

		const auto frame_count = curv_abs_frames.size();
		const auto section_count = curv_abs_section.size();

		m_indices_vectors.clear();
		m_indices_vectors.reserve(frame_count);
		m_frames_length_vectors.reserve(frame_count);
		m_beam_length_vectors.reserve(section_count);
		m_indices_vectors_draw.reserve(frame_count); // just for drawing


		/*
		 * This main loop iterates through the frames, comparing their curvilinear abscissa values with those of the
		 beam sections: If the frame's abscissa is less than the current section's, it assigns the current section
		 index. If they're equal, it assigns the current index and then increments it. If the frame's abscissa is
		 greater, it increments the index and then assigns it.
		 * */
		constexpr auto epsilon = std::numeric_limits<SReal>::epsilon();
		auto current_section_index = 1;
		for (auto i = 0; i < frame_count; ++i) {
			// The frame is associated with the current section
			if (curv_abs_section[current_section_index] > curv_abs_frames[i]) {
				m_indices_vectors.emplace_back(current_section_index);
				m_indices_vectors_draw.emplace_back(current_section_index);
			}
			// The frame is on the current section
			else if (std::abs(curv_abs_section[current_section_index] - curv_abs_frames[i]) < epsilon) {
				m_indices_vectors.emplace_back(current_section_index);
				current_section_index++;
				m_indices_vectors_draw.emplace_back(current_section_index);
			}
			// The frame is after the current section
			else {
				current_section_index++;
				m_indices_vectors.emplace_back(current_section_index);
				m_indices_vectors_draw.emplace_back(current_section_index);
			}

			// Fill the vector m_framesLengthVectors with the distance
			// between frame(output) and the closest beam node toward the base
			m_frames_length_vectors.emplace_back(curv_abs_frames[i] - curv_abs_section[m_indices_vectors.back() - 1]);
		}

		// compute the length of each beam segment.
		std::adjacent_difference(curv_abs_section.begin() + 1, curv_abs_section.end(),
								 std::back_inserter(m_beam_length_vectors));

		msg_info("BaseCosseratMapping") << "m_indicesVectors : " << m_indices_vectors << msgendl;
		std::cout << "--------------------------------------" << std::endl;
	}


	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::computeExponentialSE3(const double &sub_section_length,
																	  const Coord1 &strain_n, Frame &g_X_n) {
		const auto I4 = Mat4x4::Identity();

		// Get the angular part of the strain
		Vec3 k = Vec3(strain_n(0), strain_n(1), strain_n(2));
		SReal theta = k.norm();

		SE3 _g_X;
		SE3 Xi_hat_n = buildXiHat(strain_n);

		// todo: change double to Real
		if (theta <= std::numeric_limits<double>::epsilon()) {
			_g_X = I4 + sub_section_length * Xi_hat_n;
		} else {
			double scalar1 = (1.0 - std::cos(sub_section_length * theta)) / std::pow(theta, 2);
			double scalar2 = (sub_section_length * theta - std::sin(sub_section_length * theta)) / std::pow(theta, 3);
			// Taylor expansion of exponential
			_g_X = I4 + sub_section_length * Xi_hat_n + scalar1 * Xi_hat_n * Xi_hat_n +
				   scalar2 * Xi_hat_n * Xi_hat_n * Xi_hat_n;
		}

		Mat3x3 M;
		_g_X.getsub(0, 0, M); // get the rotation matrix

		// convert the rotation 3x3 matrix to a quaternion
		Quat<SReal> R;
		R.fromMatrix(M);
		g_X_n = Frame(Vec3(_g_X(0, 3), _g_X(1, 3), _g_X(2, 3)), R);
		// std::cout << "Translation :"<< Vec3(_g_X(0, 3), _g_X(1, 3), _g_X(2,3)) << std::endl;
		// std::cout << " ==> R : "<< R << std::endl;
	}

	// Fill exponential vectors
	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::updateExponentialSE3(const vector<Coord1> &strain_state) {
		auto curv_abs_frames = getReadAccessor(d_curv_abs_frames);


		m_frames_exponential_se3_vectors.clear();
		m_nodes_exponential_se3_vectors.clear();
		m_nodes_logarithm_se3_vectors.clear();

		const auto sz = curv_abs_frames.size();
		// Compute exponential at each frame point
		for (auto i = 0; i < sz; ++i) {
			Frame g_X_frame_i;

			const Coord1 strain_n = strain_state[m_indices_vectors[i] - 1]; // Cosserat reduce coordinates (strain)

			// the size varies from 3 to 6
			// The distance between the frame node and the closest beam node toward the base
			const SReal sub_section_length = m_frames_length_vectors[i];
			computeExponentialSE3(sub_section_length, strain_n, g_X_frame_i);
			m_frames_exponential_se3_vectors.push_back(g_X_frame_i);

			msg_info() << "_________________" << i << "_________________________" << msgendl
					   << "x :" << sub_section_length << "; strain :" << strain_n << msgendl
					   << "m_framesExponentialSE3Vectors :" << g_X_frame_i;
		}

		// Compute the exponential on the nodes
		m_nodes_exponential_se3_vectors.push_back(Frame(Vec3(0.0, 0.0, 0.0), Quat(0., 0., 0., 1.))); // The first node.
		// todo : merge this section with the previous one
		for (unsigned int j = 0; j < strain_state.size(); ++j) {
			Coord1 strain_n = strain_state[j];
			const SReal section_length = m_beam_length_vectors[j];

			Frame g_X_node_j;
			computeExponentialSE3(section_length, strain_n, g_X_node_j);
			m_nodes_exponential_se3_vectors.push_back(g_X_node_j);

			msg_info() << "_________________Beam Node Expo___________________" << msgendl
					   << "Node m_framesExponentialSE3Vectors :" << g_X_node_j << msgendl
					   << "_________________Beam Node Expo___________________";
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::compute_matrix_Adjoint(const Frame &frame, TangentTransform &adjoint) {
		Mat3x3 R = extract_rotation_matrix(frame); // extract rotation matrix frame
		Vec3 u = frame.getOrigin(); // get the linear part vec3

		Mat3x3 tilde_u_R = extract_tild_matrix(u) * R;

		adjoint.setsub(0, 0, R);
		adjoint.setsub(3, 3, R);
		adjoint.setsub(3, 0, tilde_u_R);
	}

	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::compute_matrix_CoAdjoint(const Frame &frame, Mat6x6 &coAdjoint) {
		coAdjoint.clear();
		Mat3x3 R = extract_rotation_matrix(frame); // extract rotation matrix frame
		Vec3 u = frame.getOrigin(); // get the linear part vec3
		Mat3x3 tilde_u_R = extract_tild_matrix(u) * R;

		coAdjoint.setsub(0, 0, R);
		coAdjoint.setsub(3, 3, R);
		coAdjoint.setsub(0, 3, tilde_u_R);
	}

	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::compute_matrix_adj(const Vec6 &Xi, Mat6x6 &adjoint) {
		Mat3x3 tilde_rot = extract_tild_matrix(Vec3(Xi[0], Xi[1], Xi[2]));
		Mat3x3 tilde_trans = extract_tild_matrix(Vec3(Xi[3], Xi[4], Xi[5]));

		adjoint.setsub(0, 0, tilde_rot);
		adjoint.setsub(3, 3, tilde_rot);
		adjoint.setsub(3, 0, tilde_trans);
	}

	template<class TIn1, class TIn2, class TOut>
	inline auto BaseCosseratMapping<TIn1, TIn2, TOut>::extract_tild_matrix(const Vec3 &u) -> SO3 {
		SO3 tild;

		tild[0][1] = -u[2];
		tild[0][2] = u[1];
		tild[1][2] = -u[0];
		tild[1][0] = -tild[0][1];
		tild[2][0] = -tild[0][2];
		tild[2][1] = -tild[1][2];

		return tild;
	}

	template<class TIn1, class TIn2, class TOut>
	inline auto BaseCosseratMapping<TIn1, TIn2, TOut>::buildXiHat(const Vec3 &strain_i) -> SE3 {
		SE3 xi_hat;
		xi_hat.setsub(0, 0, extract_tild_matrix(strain_i));
		// To keep the length no null,
		// This is on 0, because the beam is defined along x
		xi_hat[0][3] = 1.0;
		return xi_hat;
	}

	template<class TIn1, class TIn2, class TOut>
	inline auto BaseCosseratMapping<TIn1, TIn2, TOut>::buildXiHat(const Vec6 &strain_i) -> SE3 {
		SE3 xi_hat;
		xi_hat.setsub(0, 0, extract_tild_matrix(Vec3(strain_i(0), strain_i(1), strain_i(2))));

		for (unsigned int i = 0; i < 3; i++) {
			xi_hat[i][3] += strain_i(i + 3);
			if (xi_hat[0][3] < 0.0001)
				xi_hat[0][3] = 1.0;
		}

		return xi_hat;
	}
	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::compute_matrix_coadj(const Vec6 &Xi, Mat6x6 &adjoint) {
		Mat3x3 tilde_rot = extract_tild_matrix(Vec3(Xi[0], Xi[1], Xi[2]));
		Mat3x3 tilde_trans = extract_tild_matrix(Vec3(Xi[3], Xi[4], Xi[5]));

		adjoint.setsub(0, 0, tilde_rot);
		adjoint.setsub(3, 3, tilde_rot);
		adjoint.setsub(0, 3, tilde_trans);
	}


	template<class TIn1, class TIn2, class TOut>
	auto BaseCosseratMapping<TIn1, TIn2, TOut>::computeLogarithm(const double &x, const Mat4x4 &gX) -> Mat4x4 {
		// Compute theta before everything
		const double theta = computeTheta(x, gX);
		Mat4x4 I4 = Mat4x4::Identity();
		Mat4x4 log_gX;

		double csc_theta = 1.0 / (sin(x * theta / 2.0));
		double sec_theta = 1.0 / (cos(x * theta / 2.0));
		double cst = (1.0 / 8) * (csc_theta * csc_theta * csc_theta) * sec_theta;
		double x_theta = x * theta;
		double cos_2XTheta = cos(2.0 * x_theta);
		double cos_XTheta = cos(x_theta);
		double sin_2XTheta = sin(2.0 * x_theta);
		double sin_XTheta = sin(x_theta);

		if (theta <= std::numeric_limits<double>::epsilon())
			log_gX = I4;
		else {
			log_gX =
					cst * ((x_theta * cos_2XTheta - sin_XTheta) * I4 -
						   (x_theta * cos_XTheta + 2.0 * x_theta * cos_2XTheta - sin_XTheta - sin_2XTheta) * gX +
						   (2.0 * x_theta * cos_XTheta + x_theta * cos_2XTheta - sin_XTheta - sin_2XTheta) * (gX * gX) -
						   (x_theta * cos_XTheta - sin_XTheta) * (gX * gX * gX));
		}

		return log_gX;
	}

	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::updateTangExpSE3(const vector<Coord1> &inDeform) {

		// Curv abscissa of nodes and frames
		auto curv_abs_section = getReadAccessor(d_curv_abs_section);
		auto curv_abs_frames = getReadAccessor(d_curv_abs_frames);

		unsigned int sz = curv_abs_frames.size();
		m_frames_tang_exp_vectors.resize(sz);

		// Compute tangExpo at frame points
		for (unsigned int i = 0; i < sz; i++) {
			TangentTransform tangent_gX_frame;

			Coord1 strain_frame_i = inDeform[m_indices_vectors[i] - 1];
			double curv_abs_x_i = m_frames_length_vectors[i];
			//computeTangExp(curv_abs_x_i, strain_frame_i, tangent_gX_frame);

			if constexpr (Coord1::static_size == 3)
				computeTangExpImplementation(curv_abs_x_i,
											 Vec6(strain_frame_i(0), strain_frame_i(1), strain_frame_i(2), 0, 0, 0),
											 tangent_gX_frame);
			else
				computeTangExpImplementation(curv_abs_x_i, strain_frame_i, tangent_gX_frame);

			m_frames_tang_exp_vectors[i] = tangent_gX_frame;

			msg_info() << "x :" << curv_abs_x_i << "; k :" << strain_frame_i << msgendl
					   << "m_framesTangExpVectors :" << m_frames_tang_exp_vectors[i];
		}

		// Compute the TangExpSE3 at the nodes
		m_nodes_tang_exp_vectors.clear();
		TangentTransform tangExpO;
		tangExpO.clear();
		m_nodes_tang_exp_vectors.push_back(tangExpO);

		for (size_t j = 1; j < curv_abs_section.size(); j++) {
			Coord1 strain_node_i = inDeform[j - 1];
			double curv_abs_node_x = m_beam_length_vectors[j - 1];
			TangentTransform tangent_gX_node;
			tangent_gX_node.clear();

			if constexpr (Coord1::static_size == 3)
				computeTangExpImplementation(curv_abs_node_x,
											 Vec6(strain_node_i(0), strain_node_i(1), strain_node_i(2), 0, 0, 0),
											 tangent_gX_node);
			else
				computeTangExpImplementation(curv_abs_node_x, strain_node_i, tangent_gX_node);

			m_nodes_tang_exp_vectors.push_back(tangent_gX_node);
		}
		msg_info() << "Node TangExpo : " << m_nodes_tang_exp_vectors;
	}

	// template <class TIn1, class TIn2, class TOut>
	// void BaseCosseratMapping<TIn1, TIn2, TOut>::computeTangExp(double &curv_abs_n,
	//                                                            const Coord1 &strain_i,
	//                                                            Mat6x6 &TgX)
	// {
	//     if constexpr( Coord1::static_size == 3 )
	//         computeTangExpImplementation(curv_abs_n, Vec6(strain_i(0),strain_i(1),strain_i(2),0,0,0), TgX);
	//     else
	//         computeTangExpImplementation(curv_abs_n, strain_i, TgX);
	// }

	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::computeTangExpImplementation(double &curv_abs_n, const Vec6 &strain_i,
																			 Mat6x6 &TgX) {
		SReal theta = Vec3(strain_i(0), strain_i(1), strain_i(2)).norm();
		SO3 tilde_k = extract_tild_matrix(Vec3(strain_i(0), strain_i(1), strain_i(2)));
		SO3 tilde_q = extract_tild_matrix(Vec3(strain_i(3), strain_i(4), strain_i(5)));

		Mat6x6 ad_Xi;
		ad_Xi.setsub(0, 0, tilde_k);
		ad_Xi.setsub(3, 3, tilde_k);
		ad_Xi.setsub(3, 0, tilde_q);


		Mat6x6 Id6 = Mat6x6::Identity();
		if (theta <= std::numeric_limits<double>::epsilon()) {
			double scalar0 = std::pow(curv_abs_n, 2) / 2.0;
			TgX = curv_abs_n * Id6 + scalar0 * ad_Xi;
		} else {
			double scalar1 = (4.0 - 4.0 * cos(curv_abs_n * theta) - curv_abs_n * theta * sin(curv_abs_n * theta)) /
							 (2.0 * theta * theta);
			double scalar2 = (4.0 * curv_abs_n * theta + curv_abs_n * theta * cos(curv_abs_n * theta) -
							  5.0 * sin(curv_abs_n * theta)) /
							 (2.0 * theta * theta * theta);
			double scalar3 = (2.0 - 2.0 * cos(curv_abs_n * theta) - curv_abs_n * theta * sin(curv_abs_n * theta)) /
							 (2.0 * theta * theta * theta * theta);
			double scalar4 = (2.0 * curv_abs_n * theta + curv_abs_n * theta * cos(curv_abs_n * theta) -
							  3.0 * sin(curv_abs_n * theta)) /
							 (2.0 * theta * theta * theta * theta * theta);

			TgX = curv_abs_n * Id6 + scalar1 * ad_Xi + scalar2 * ad_Xi * ad_Xi + scalar3 * ad_Xi * ad_Xi * ad_Xi +
				  scalar4 * ad_Xi * ad_Xi * ad_Xi * ad_Xi;
		}
	}

	// template<class TIn1, class TIn2, class TOut>
	// [[maybe_unused]] Vec6 BaseCosseratMapping<TIn1, TIn2, TOut>::computeETA(const Vec6 &baseEta,
	// 																		const vector<Deriv1> &k_dot,
	// 																		const double abs_input) {
	// 	// Get the positions from model 0. This function returns the position wrapped in a Data<>
	// 	auto d_x1 = m_strain_state->read(sofa::core::vec_id::read_access::position);
	//
	// 	// To access the actual content (in this case position) from a data, we have to use
	// 	// a read accessor that insures the data is updated according to DDGNode state
	// 	auto x1 = getReadAccessor(*d_x1);
	//
	// 	// Same as for x1, query a read accessor so we can access the content of d_curv_abs_section
	// 	auto curv_abs_input = getReadAccessor(d_curv_abs_section);
	//
	// 	auto &kdot = k_dot[m_index_input];
	// 	Vec6 Xi_dot{kdot[0], kdot[1], kdot[2], 0, 0, 0};
	//
	// 	// if m_indexInput is == 0
	// 	double diff0 = abs_input;
	// 	double _diff0 = -abs_input;
	//
	// 	if (m_index_input != 0) {
	// 		diff0 = abs_input - curv_abs_input[m_index_input - 1];
	// 		_diff0 = curv_abs_input[m_index_input - 1] - abs_input;
	// 	}
	//
	// 	Frame outTransform;
	// 	computeExponentialSE3(_diff0, x1[m_index_input], outTransform);
	//
	// 	TangentTransform adjointMatrix;
	// 	compute_matrix_Adjoint(outTransform, adjointMatrix);
	//
	// 	TangentTransform tangentMatrix;
	//
	// 	computeTangExp(diff0, x1[m_index_input], tangentMatrix);
	//
	// 	return adjointMatrix * (baseEta + tangentMatrix * Xi_dot);
	// }


	template<class TIn1, class TIn2, class TOut>
	double BaseCosseratMapping<TIn1, TIn2, TOut>::computeTheta(const double &x, const Mat4x4 &gX) {
		double Tr_gx = sofa::type::trace(gX);

		if (x > std::numeric_limits<double>::epsilon())
			return (1.0 / x) * std::acos((Tr_gx / 2.0) - 1);

		return 0.0;
	}

	template<class TIn1, class TIn2, class TOut>
	void BaseCosseratMapping<TIn1, TIn2, TOut>::printMatrix(const Mat6x6 R) {
		// TODO(dmarchal: 2024/06/07): Remove the use of printf in addition to
		// reconsider the implementation of common utility functions in instance
		// method.
		for (unsigned int k = 0; k < 6; k++) {
			for (unsigned int i = 0; i < 6; i++)
				printf("  %lf", R[k][i]);
			printf("\n");
		}
	}

	template<class TIn1, class TIn2, class TOut>
	Mat3x3 BaseCosseratMapping<TIn1, TIn2, TOut>::extract_rotation_matrix(const Frame &frame) {

		Quat q = frame.getOrientation();
		Mat3x3 mat33;
		q.toMatrix(mat33);
		return mat33;
	}

	template<class TIn1, class TIn2, class TOut>
	auto BaseCosseratMapping<TIn1, TIn2, TOut>::buildProjector(const Frame &T) -> TangentTransform {
		TangentTransform P; // It's a 6x6 matrix

		Mat3x3 mat33 = T.getRotationMatrix();

		P.setsub(0,3,mat33);
		P.setsub(3,0, mat33);
		return P;
	}


	template<class TIn1, class TIn2, class TOut>
	auto BaseCosseratMapping<TIn1, TIn2, TOut>::convertTransformToMatrix4x4(const Frame &T) -> Mat4x4 {
		Mat4x4 M = Mat4x4::Identity();
		Mat3x3 R = extract_rotation_matrix(T); // extract the rotation matrix from frame (quat,vec3)
		Vec3 trans = T.getOrigin();
		M.setsub(0, 0, R);
		M.setsub(0, 3, trans);

		return M;
	}

	template<class TIn1, class TIn2, class TOut>
	auto BaseCosseratMapping<TIn1, TIn2, TOut>::piecewiseLogmap(const _SE3 &g_x) -> Vec6 {
		_SE3 Xi_hat;

		double x = 1.0;
		double theta = std::acos(g_x.trace() / 2.0 - 1.0);

		if (theta == 0) {
			Xi_hat = 1.0 / x * (g_x - Matrix4d::Identity());
		} else {
			double x_theta = x * theta;
			double sin_x_theta = std::sin(x_theta);
			double cos_x_theta = std::cos(x_theta);
			double t3 = 2 * sin_x_theta * cos_x_theta;
			double t4 = 1 - 2 * sin_x_theta * sin_x_theta;
			double t5 = x_theta * t4;

			Matrix4d gp2 = g_x * g_x;
			Matrix4d gp3 = gp2 * g_x;

			Xi_hat = 1.0 / x *
					 (0.125 * (1.0 / std::sin(x_theta / 2.0) / std::sin(x_theta / 2.0) / std::sin(x_theta / 2.0)) *
					  std::cos(x_theta / 2.0) *
					  ((t5 - sin_x_theta) * Matrix4d::Identity() -
					   (x_theta * cos_x_theta + 2 * t5 - sin_x_theta - t3) * g_x +
					   (2 * x_theta * cos_x_theta + t5 - sin_x_theta - t3) * gp2 -
					   (x_theta * cos_x_theta - sin_x_theta) * gp3));
		}

		Vec6 xci = Vec6(Xi_hat(2, 1), Xi_hat(0, 2), Xi_hat(1, 0), Xi_hat(0, 3), Xi_hat(1, 3), Xi_hat(2, 3));
		return xci;
	}

} // namespace Cosserat::mapping

