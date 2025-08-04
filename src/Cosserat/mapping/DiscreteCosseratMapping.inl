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
#include <Cosserat/mapping/DiscreteCosseratMapping.h>

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/gl/template.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/helper/visual/DrawTool.h>
#include <sofa/type/Quat.h>

#include <string>

namespace Cosserat::mapping {

	using sofa::core::objectmodel::BaseContext;
	using sofa::defaulttype::SolidTypes;
	using sofa::helper::AdvancedTimer;
	using sofa::helper::WriteAccessor;
	using sofa::type::RGBAColor;

	template<class TIn1, class TIn2, class TOut>
	DiscreteCosseratMapping<TIn1, TIn2, TOut>::DiscreteCosseratMapping() :
		d_deformationAxis(initData(&d_deformationAxis, (int) 1, "deformationAxis",
								   "the axis in which we want to show the deformation.\n")),
		d_max(initData(&d_max, (SReal) 1.0e-2, "max", "the maximum of the deformation.\n")),
		d_min(initData(&d_min, (SReal) 0.0, "min", "the minimum of the deformation.\n")),
		d_radius(initData(&d_radius, (SReal) 0.05, "radius", "the axis in which we want to show the deformation.\n")),
		d_drawMapBeam(initData(&d_drawMapBeam, true, "nonColored",
							   "if this parameter is false, you draw the beam with "
							   "color according to the force apply to each beam")),
		d_color(initData(&d_color, sofa::type::RGBAColor(40 / 255.0, 104 / 255.0, 137 / 255.0, 0.8), "color",
						 "The default beam color")),
		d_index(initData(&d_index, "index",
						 "if this parameter is false, you draw the beam with color "
						 "according to the force apply to each beam")),
		d_baseIndex(initData(&d_baseIndex, (unsigned int) 0, "baseIndex",
							 "This parameter defines the index of the rigid "
							 "base of Cosserat models, 0 by default this can"
							 "take another value if the rigid base is given "
							 "by another body.")) {
		this->addUpdateCallback("updateFrames", {&d_curv_abs_section, &d_curv_abs_frames, &d_debug},
								[this](const sofa::core::DataTracker &t) {
									SOFA_UNUSED(t);
									this->update_geometry_info();

									const sofa::VecCoord_t<In1> &strain_state =
											m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

									this->updateExponentialSE3(strain_state);
									return sofa::core::objectmodel::ComponentState::Valid;
								},
								{});
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::doBaseCosseratInit() {
		m_colorMap.setColorScheme("Blue to Red");
		m_colorMap.reinit();
	}

	template<class TIn1, class TIn2, class TOut>
	void
	DiscreteCosseratMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams * /* mparams */,
													 const vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
													 const vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
													 const vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) {

		if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		// Checking the componentState, to trigger a callback if other data fields (specifically
		// d_curv_abs_section and d_curv_abs_frames) were changed dynamically
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;
		/// Do Apply
		// We need only one input In model and input Root model (if present)
		const sofa::VecCoord_t<In1> &in1 = dataVecIn1Pos[0]->getValue();
		const sofa::VecCoord_t<In2> &in2 = dataVecIn2Pos[0]->getValue();

		const auto sz = d_curv_abs_frames.getValue().size();
		sofa::VecCoord_t<Out> &out = *dataVecOutPos[0]->beginEdit(); // frames states
		out.resize(sz);
		const auto baseIndex = d_baseIndex.getValue();

		// update the Exponential matrices according to new deformation
		// Here we update m_framesExponentialSE3Vectors & m_nodesExponentialSE3Vectors
		// Which are the homogeneous matrices of the frames and the nodes in local
		// coordinate.
		this->updateExponentialSE3(in1);

		/* Apply the transformation to go from cossserat to SOFA frame*/
		const auto frame0 = Frame(In2::getCPos(in2[baseIndex]), In2::getCRot(in2[baseIndex]));

		// Cache the printLog value out of the loop, otherwise it will trigger a graph
		// update at every iteration.
		bool doPrintLog = this->f_printLog.getValue();
		for (unsigned int i = 0; i < sz; i++) {
			auto frame = frame0;
			std::cout << "frame "<< i <<std::endl;
			for (unsigned int u = 0; u < m_indices_vectors[i]; u++) {
				std::cout << "section "<< u << std::endl;
				std::cout << "Expo node : " << m_nodes_exponential_se3_vectors[u] <<std::endl;
				std::cout<< "current frame : "<< frame << std::endl;
				frame *= m_nodes_exponential_se3_vectors[u]; // frame = gX(L_0)*...*gX(L_{n-1})
			}
			frame *= m_frames_exponential_se3_vectors[i]; // frame*gX(x)

			// This is a lazy printing approach, so there is no time consuming action in
			// the core of the loop.
			msg_info_when(doPrintLog) << "Frame  : " << i << " = " << frame;

			Vec3 origin = frame.getOrigin();
			Quat orientation = frame.getOrientation();
			out[i] = sofa::Coord_t<Out>(origin, orientation);
		}

		// If the printLog attribute is checked then print distance between out
		// frames.
		if (doPrintLog) {
			std::stringstream tmp;
			for (unsigned int i = 0; i < out.size() - 1; i++) {
				Vec3 diff = out[i + 1].getCenter() - out[i].getCenter();
				tmp << "dist " << i << "  : " << diff.norm() << msgendl;
			}
			msg_info() << tmp.str();
		}

		// Debug output if needed
		if (this->f_printLog.getValue()) {
			displayOutputFrames(out, "apply - computed output frames");
			displayTransformMatrices("apply - transformation matrices");
		}

		// TODO(dmarchal:2024/06/13): This looks a suspicious design pattern,
		// elaborate more on the purpose of m_indexInput and how to use it.
		m_index_input = 0;
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::computeLogarithm(const double &x, const Mat4x4 &gX,
																	 Mat4x4 &log_gX) {

		// computes theta
		const double theta = computeTheta(x, gX);

		// if theta is very small, we return log_gX as the identity.
		if (theta <= std::numeric_limits<double>::epsilon()) {
			log_gX = Mat4x4::Identity();
			return;
		}

		// otherwise we compute it
		const double csc_theta = 1.0 / (sin(x * theta / 2.0));
		const double sec_theta = 1.0 / (cos(x * theta / 2.0));
		const double cst = (1.0 / 8) * (csc_theta * csc_theta * csc_theta) * sec_theta;
		const double x_theta = x * theta;
		const double cos_2x_theta = cos(2.0 * x_theta);
		const double cos_x_theta = cos(x_theta);
		const double sin_2x_theta = sin(2.0 * x_theta);
		const double sin_x_theta = sin(x_theta);

		log_gX.clear();
		log_gX =
				cst * ((x_theta * cos_2x_theta - sin_x_theta) * Mat4x4::Identity() -
					   (x_theta * cos_x_theta + 2.0 * x_theta * cos_2x_theta - sin_x_theta - sin_2x_theta) * gX +
					   (2.0 * x_theta * cos_x_theta + x_theta * cos_2x_theta - sin_x_theta - sin_2x_theta) * (gX * gX) -
					   (x_theta * cos_x_theta - sin_x_theta) * (gX * gX * gX));
	}

	template<class TIn1, class TIn2, class TOut>
	void
	DiscreteCosseratMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams * /* mparams */,
													  const vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
													  const vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
													  const vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) {

		if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;
		if (d_debug.getValue())
			std::cout << " ########## ApplyJ Function ########" << std::endl;
		const sofa::VecDeriv_t<In1> &in1_vel = dataVecIn1Vel[0]->getValue();
		const sofa::VecDeriv_t<In2> &in2_vel = dataVecIn2Vel[0]->getValue();
		sofa::VecDeriv_t<Out> &out_vel = *dataVecOutVel[0]->beginEdit();
		const auto baseIndex = d_baseIndex.getValue();

		// Curv abscissa of nodes and frames
		sofa::helper::ReadAccessor<sofa::Data<vector<double>>> curv_abs_section = d_curv_abs_section;
		sofa::helper::ReadAccessor<sofa::Data<vector<double>>> curv_abs_frames = d_curv_abs_frames;

		const sofa::VecDeriv_t<In1> &inDeform =
				m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue(); // strains
		//1. Compute the tangent Exponential SE3 vectors
		this->updateTangExpSE3(inDeform);

		//2. Get base velocity and convert to Vec6, for the facility of computation
		// Get base velocity as input this is also called eta
		//2.1 Get the base velocity from input
		Vec6 baseVelocity; //
		for (auto u = 0; u < 6; u++)
			baseVelocity[u] = in2_vel[baseIndex][u];

		//2.2 Apply the local transform i.e. from SOFA's frame to Cosserat's frame
		const sofa::VecCoord_t<In2> &xfrom2Data =
				m_rigid_base->read(sofa::core::vec_id::read_access::position)->getValue();
		// Get the transformation from the SOFA to the local frame
		auto TInverse = Frame(xfrom2Data[baseIndex].getCenter(), xfrom2Data[baseIndex].getOrientation()).inversed();
		Mat6x6 P = this->buildProjector(TInverse);

		// List of velocity vectors at nodes
		m_nodes_velocity_vectors.clear();
		Vec6 baseLocalVelocity = P * baseVelocity; // This is the base velocity in Locale frame

		m_nodes_velocity_vectors.push_back(baseLocalVelocity);
		if (d_debug.getValue())
			std::cout << "Base local Velocity :" << baseLocalVelocity << std::endl;

		// Compute velocity at nodes
		for (unsigned int i = 1; i < curv_abs_section.size(); i++) {
			auto Trans = m_nodes_exponential_se3_vectors[i].inversed();
			TangentTransform Adjoint;
			Adjoint.clear();
			this->computeAdjoint(Trans, Adjoint);

			/// The null vector is replace by the linear velocity in Vec6Type
			Vec6 Xi_dot = Vec6(in1_vel[i - 1], Vec3(0.0, 0.0, 0.0));

			Vec6 eta_node_i = Adjoint * (m_nodes_velocity_vectors[i - 1] + m_nodes_tang_exp_vectors[i] * Xi_dot);
			m_nodes_velocity_vectors.push_back(eta_node_i);
			if (d_debug.getValue())
				std::cout << "Node velocity : " << i << " = " << eta_node_i << std::endl;
		}

		const sofa::VecCoord_t<Out> &out = m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		auto sz = curv_abs_frames.size();
		out_vel.resize(sz);
		for (unsigned int i = 0; i < sz; i++) {
			auto Trans = m_frames_exponential_se3_vectors[i].inversed();
			TangentTransform Adjoint; ///< the class insure that the constructed adjoint is zeroed.
			Adjoint.clear();
			this->computeAdjoint(Trans, Adjoint);
			Vec6 frame_Xi_dot;

			for (auto u = 0; u < 3; u++) {
				frame_Xi_dot(u) = in1_vel[m_indices_vectors[i] - 1][u];
				frame_Xi_dot(u + 3) = 0.;
			}
			Vec6 eta_frame_i = Adjoint * (m_nodes_velocity_vectors[m_indices_vectors[i] - 1] +
										  m_frames_tang_exp_vectors[i] * frame_Xi_dot); // eta

			auto T = Frame(out[i].getCenter(), out[i].getOrientation());
			Mat6x6 Proj = this->buildProjector(T);

			out_vel[i] = Proj * eta_frame_i;

			if (d_debug.getValue())
				std::cout << "Frame velocity : " << i << " = " << eta_frame_i << std::endl;
		}
		
		// Debug output if needed
		if (this->f_printLog.getValue()) {
			displayInputVelocities(in1_vel, in2_vel, "applyJ - input velocities");
			displayOutputVelocities(out_vel, "applyJ - computed output velocities");
		}
		
		dataVecOutVel[0]->endEdit();
		m_index_input = 0;
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			const vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
			const vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
			const vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) {

		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## ApplyJT force Function ########" << std::endl;
		const sofa::VecDeriv_t<Out> &in = dataVecInForce[0]->getValue();

		sofa::VecDeriv_t<In1> &out1 = *dataVecOut1Force[0]->beginEdit();
		sofa::VecDeriv_t<In2> &out2 = *dataVecOut2Force[0]->beginEdit();
		const auto baseIndex = d_baseIndex.getValue();

		const sofa::VecCoord_t<Out> &frame =
				m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		const sofa::DataVecCoord_t<In1> *x1fromData = m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> x1from = x1fromData->getValue();
		vector<Vec6> local_F_Vec;
		local_F_Vec.clear();

		out1.resize(x1from.size());

		// convert the input from Deriv type to vec6 type, for the purpose of the
		// matrix vector multiplication
		for (unsigned int var = 0; var < in.size(); ++var) {
			Vec6 vec;
			for (unsigned j = 0; j < 6; j++)
				vec[j] = in[var][j];
			// Convert input from global frame(SOFA) to local frame
			const auto _T = Frame(frame[var].getCenter(), frame[var].getOrientation());
			Mat6x6 P_trans = (this->buildProjector(_T));
			P_trans.transpose();
			Vec6 local_F = P_trans * vec;
			local_F_Vec.push_back(local_F);
		}

		// Compute output forces
		auto sz = m_indices_vectors.size();
		auto index = m_indices_vectors[sz - 1];
		m_total_beam_force_vectors.clear();
		m_total_beam_force_vectors.resize(sz);

		Vec6 F_tot;
		F_tot.clear();
		m_total_beam_force_vectors.push_back(F_tot);

		Mat3x6 matB_trans;
		matB_trans.clear();
		for (unsigned int k = 0; k < 3; k++)
			matB_trans[k][k] = 1.0;

		for (auto s = sz; s--;) {
			Mat6x6 coAdjoint;

			this->computeCoAdjoint(m_frames_exponential_se3_vectors[s],
								   coAdjoint); // m_framesExponentialSE3Vectors[s] computed in apply
			Vec6 node_F_Vec = coAdjoint * local_F_Vec[s];
			Mat6x6 temp = m_frames_tang_exp_vectors[s]; // m_framesTangExpVectors[s] computed in
			// applyJ (here we transpose)
			temp.transpose();
			Vec3 f = matB_trans * temp * node_F_Vec;

			if (index != m_indices_vectors[s]) {
				index--;
				// bring F_tot to the reference of the new beam
				this->computeCoAdjoint(m_nodes_exponential_se3_vectors[index],
									   coAdjoint); // m_nodes_exponential_se3_vectors computed in apply
				F_tot = coAdjoint * F_tot;
				Mat6x6 temp = m_nodes_tang_exp_vectors[index];
				temp.transpose();
				// apply F_tot to the new beam
				Vec3 temp_f = matB_trans * temp * F_tot;
				out1[index - 1] += temp_f;
			}
			if (d_debug.getValue())
				std::cout << "f at s =" << s << " and index" << index << " is : " << f << std::endl;

			// compute F_tot
			F_tot += node_F_Vec;
			out1[m_indices_vectors[s] - 1] += f;
		}

		auto frame0 = Frame(frame[0].getCenter(), frame[0].getOrientation());
		Mat6x6 M = this->buildProjector(frame0);
		out2[baseIndex] += M * F_tot;

		if (d_debug.getValue()) {
			std::cout << "Node forces " << out1 << std::endl;
			std::cout << "base Force: " << out2[baseIndex] << std::endl;
		}

		// Debug output if needed
		if (this->f_printLog.getValue()) {
			displayOutputForces(in, "applyJT - input forces");
			displayInputForces(out1, out2, "applyJT - computed input forces");
		}

		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::ConstraintParams * /*cparams*/,
			const vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
			const vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
			const vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) {
		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## ApplyJT Constraint Function ########" << std::endl;
		// We need only one input In model and input Root model (if present)
		sofa::MatrixDeriv_t<In1> &out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the strain space
		// (reduced coordinate)
		// constraints on the reference frame (base frame)
		sofa::MatrixDeriv_t<In2> &out2 = *dataMatOut2Const[0]->beginEdit();
		// constraints on the reference frame (base frame)
		const sofa::MatrixDeriv_t<Out> &in = dataMatInConst[0]->getValue();

		const sofa::VecCoord_t<Out> &frame =
				m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		const sofa::DataVecCoord_t<In1> *x1fromData = m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> x1from = x1fromData->getValue();

		Mat3x6 matB_trans;
		matB_trans.clear();
		for (unsigned int k = 0; k < 3; k++)
			matB_trans[k][k] = 1.0;

		vector<std::tuple<int, Vec6>> NodesInvolved;
		vector<std::tuple<int, Vec6>> NodesInvolvedCompressed;
		// helper::vector<Vec6> NodesConstraintDirection;

		typename sofa::MatrixDeriv_t<Out>::RowConstIterator rowItEnd = in.end();

		for (typename sofa::MatrixDeriv_t<Out>::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt) {
			if (d_debug.getValue()) {
				std::cout << "************* Apply JT (MatrixDeriv) iteration on line ";
				std::cout << rowIt.index();
				std::cout << "*************  " << std::endl;
			}
			typename sofa::MatrixDeriv_t<Out>::ColConstIterator colIt = rowIt.begin();
			typename sofa::MatrixDeriv_t<Out>::ColConstIterator colItEnd = rowIt.end();

			// Creates a constraints if the input constraint is not empty.
			if (colIt == colItEnd) {
				if (d_debug.getValue()) {
					std::cout << "no column for this constraint" << std::endl;
				}
				continue;
			}
			typename sofa::MatrixDeriv_t<In1>::RowIterator o1 =
					out1.writeLine(rowIt.index()); // we store the constraint number
			typename sofa::MatrixDeriv_t<In2>::RowIterator o2 = out2.writeLine(rowIt.index());

			NodesInvolved.clear();
			while (colIt != colItEnd) {
				int childIndex = colIt.index();

				const sofa::Deriv_t<Out> valueConst_ = colIt.val();
				Vec6 valueConst;
				for (unsigned j = 0; j < 6; j++)
					valueConst[j] = valueConst_[j];

				int indexBeam = m_indices_vectors[childIndex];

				const auto _T = Frame(frame[childIndex].getCenter(), frame[childIndex].getOrientation());
				Mat6x6 P_trans = (this->buildProjector(_T));
				P_trans.transpose();

				Mat6x6 co_adjoint;
				this->computeCoAdjoint(m_frames_exponential_se3_vectors[childIndex],
									   co_adjoint); // m_frames_exponential_se3_vectors[s] computed in apply
				Mat6x6 temp = m_frames_tang_exp_vectors[childIndex]; // m_framesTangExpVectors[s]
				// computed in applyJ
				// (here we transpose)
				temp.transpose();

				Vec6 local_F = co_adjoint * P_trans * valueConst; // constraint direction in local frame of the beam.

				Vec3 f = matB_trans * temp * local_F; // constraint direction in the strain space.

				o1.addCol(indexBeam - 1, f);
				std::tuple<int, Vec6> test = std::make_tuple(indexBeam, local_F);

				NodesInvolved.push_back(test);
				colIt++;
			}
			if (d_debug.getValue()) {
				std::cout << "==> NodesInvolved : " << std::endl;
				for (size_t i = 0; i < NodesInvolved.size(); i++)
					std::cout << "index :" << get<0>(NodesInvolved[i]) << " force :" << get<1>(NodesInvolved[i])
							  << "\n ";
			}

			// sort the Nodes Invoved by decreasing order
			std::sort(begin(NodesInvolved), end(NodesInvolved),
					  [](std::tuple<int, Vec6> const &t1, std::tuple<int, Vec6> const &t2) {
						  return std::get<0>(t1) > std::get<0>(t2); // custom compare function
					  });

			NodesInvolvedCompressed.clear();

			for (unsigned n = 0; n < NodesInvolved.size(); n++) {
				std::tuple<int, Vec6> test_i = NodesInvolved[n];
				int numNode_i = std::get<0>(test_i);
				Vec6 cumulativeF = std::get<1>(test_i);

				if (n < NodesInvolved.size() - 1) {
					std::tuple<int, Vec6> test_i1 = NodesInvolved[n + 1];
					int numNode_i1 = std::get<0>(test_i1);

					while (numNode_i == numNode_i1) {
						cumulativeF += std::get<1>(test_i1);
						//// This was if ((n!=NodesInvolved.size()-2)||(n==0)) before and I
						/// change it to
						/// if ((n!=NodesInvolved.size()-1)||(n==0)) since the code can't
						/// leave the will loop
						if ((n != NodesInvolved.size() - 1) || (n == 0)) {
							n++;
							break;
						}
						test_i1 = NodesInvolved[n + 1];
						numNode_i1 = std::get<0>(test_i1);
					}
				}
				NodesInvolvedCompressed.push_back(std::make_tuple(numNode_i, cumulativeF));
			}

			if (d_debug.getValue()) {
				std::cout << " NodesInvolved after sort and compress" << std::endl;
				for (size_t i = 0; i < NodesInvolvedCompressed.size(); i++)
					std::cout << "index :" << get<0>(NodesInvolvedCompressed[i])
							  << " force :" << get<1>(NodesInvolvedCompressed[i]) << "\n ";
			}

			for (unsigned n = 0; n < NodesInvolvedCompressed.size(); n++) {
				std::tuple<int, Vec6> test = NodesInvolvedCompressed[n];
				int numNode = std::get<0>(test);
				int i = numNode;
				Vec6 CumulativeF = std::get<1>(test);

				while (i > 0) {
					// cumulate on beam frame
					Mat6x6 coAdjoint;
					this->computeCoAdjoint(m_nodes_exponential_se3_vectors[i - 1],
										   coAdjoint); // m_nodes_exponential_se3_vectors computed in apply
					CumulativeF = coAdjoint * CumulativeF;
					// transfer to strain space (local coordinates)
					Mat6x6 temp = m_nodes_tang_exp_vectors[i - 1];
					temp.transpose();
					Vec3 temp_f = matB_trans * temp * CumulativeF;

					if (i > 1)
						o1.addCol(i - 2, temp_f);
					i--;
				}
				const auto frame0 = Frame(frame[0].getCenter(), frame[0].getOrientation());
				const Mat6x6 M = this->buildProjector(frame0);

				const Vec6 base_force = M * CumulativeF;
				o2.addCol(d_baseIndex.getValue(), base_force);
			}
		}

		//"""END ARTICULATION SYSTEM MAPPING"""
		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::computeBBox(const sofa::core::ExecParams *, bool) {
		const sofa::VecCoord_t<Out> &x = m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();

		SReal minBBox[3] = {std::numeric_limits<SReal>::max(), std::numeric_limits<SReal>::max(),
							std::numeric_limits<SReal>::max()};
		SReal maxBBox[3] = {-std::numeric_limits<SReal>::max(), -std::numeric_limits<SReal>::max(),
							-std::numeric_limits<SReal>::max()};
		for (std::size_t i = 0; i < x.size(); i++) {
			const sofa::Coord_t<Out> &p = x[i];
			for (int c = 0; c < 3; c++) {
				if (p[c] > maxBBox[c])
					maxBBox[c] = p[c];
				if (p[c] < minBBox[c])
					minBBox[c] = p[c];
			}
		}
		this->f_bbox.setValue(sofa::type::TBoundingBox<SReal>(minBBox, maxBBox));
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		if (!vparams->displayFlags().getShowMechanicalMappings())
			return;

		// draw cable
		typedef RGBAColor RGBAColor;

		const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

		const sofa::DataVecCoord_t<Out> *xfromData = m_global_frames->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<Out> xData = xfromData->getValue();
		vector<Vec3> positions;
		vector<sofa::type::Quat<SReal>> Orientation;
		positions.clear();
		Orientation.clear();
		unsigned int sz = xData.size();
		for (unsigned int i = 0; i < sz; i++) {
			positions.push_back(xData[i].getCenter());
			Orientation.push_back(xData[i].getOrientation());
		}

		// Get access articulated
		const sofa::DataVecCoord_t<In1> *artiData = m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> xPos = artiData->getValue();

		RGBAColor drawColor = d_color.getValue();
		// draw each segment of the beam as a cylinder.
		for (unsigned int i = 0; i < sz - 1; i++)
			vparams->drawTool()->drawCylinder(positions[i], positions[i + 1], d_radius.getValue(), drawColor);

		// Define color map
		SReal min = d_min.getValue();
		SReal max = d_max.getValue();
		sofa::helper::ColorMap::evaluator<SReal> _eval = m_colorMap.getEvaluator(min, max);

		glLineWidth(d_radius.getValue());
		glBegin(GL_LINES);
		if (d_drawMapBeam.getValue()) {
			sofa::type::RGBAColor _color = d_color.getValue();
			RGBAColor colorL = RGBAColor(_color[0], _color[1], _color[2], _color[3]);
			glColor4f(colorL[0], colorL[1], colorL[2], colorL[3]);
			for (unsigned int i = 0; i < sz - 1; i++) {
				vparams->drawTool()->drawLine(positions[i], positions[i + 1], colorL);
			}
		} else {
			int j = 0;
			vector<int> index = d_index.getValue();
			for (unsigned int i = 0; i < sz - 1; i++) {
				j = m_indices_vectors_draw[i] - 1; // to get the articulation on which the frame is related to
				RGBAColor color = _eval(xPos[j][d_deformationAxis.getValue()]);
				vparams->drawTool()->drawLine(positions[i], positions[i + 1], color);
			}
		}
	glLineWidth(1);
	if (!vparams->displayFlags().getShowMappings())
		if (!d_debug.getValue())
			return;
	
	// Debug output if needed
	if (this->f_printLog.getValue()) {
		displayOutputFrames(xData, "draw - rendering frames");
	}
	
	glEnd();
}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::displayOutputFrames(const sofa::VecCoord_t<Out> &frames, const std::string &label) {
		std::cout << label << std::endl;
		for (size_t i = 0; i < frames.size(); ++i) {
			std::cout << "Frame " << i << ": position=" << frames[i].getCenter() 
					  << ", orientation=" << frames[i].getOrientation() << std::endl;
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::displayInputVelocities(const sofa::VecDeriv_t<In1> &in1Vel, const sofa::VecDeriv_t<In2> &in2Vel, const std::string &label) {
		std::cout << label << std::endl;
		std::cout << "Input 1 velocities:" << std::endl;
		for (size_t i = 0; i < in1Vel.size(); ++i) {
			std::cout << "  Vel1[" << i << "]: " << in1Vel[i] << std::endl;
		}
		std::cout << "Input 2 velocities:" << std::endl;
		for (size_t i = 0; i < in2Vel.size(); ++i) {
			std::cout << "  Vel2[" << i << "]: " << in2Vel[i] << std::endl;
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::displayOutputVelocities(const sofa::VecDeriv_t<Out> &outVel, const std::string &label) {
		std::cout << label << std::endl;
		for (size_t i = 0; i < outVel.size(); ++i) {
			std::cout << "Output velocity[" << i << "]: " << outVel[i] << std::endl;
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::displayInputForces(const sofa::VecDeriv_t<In1> &in1Force, const sofa::VecDeriv_t<In2> &in2Force, const std::string &label) {
		std::cout << label << std::endl;
		std::cout << "Input 1 forces:" << std::endl;
		for (size_t i = 0; i < in1Force.size(); ++i) {
			std::cout << "  Force1[" << i << "]: " << in1Force[i] << std::endl;
		}
		std::cout << "Input 2 forces:" << std::endl;
		for (size_t i = 0; i < in2Force.size(); ++i) {
			std::cout << "  Force2[" << i << "]: " << in2Force[i] << std::endl;
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::displayOutputForces(const sofa::VecDeriv_t<Out> &outForce, const std::string &label) {
		std::cout << label << std::endl;
		for (size_t i = 0; i < outForce.size(); ++i) {
			std::cout << "Output force[" << i << "]: " << outForce[i] << std::endl;
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteCosseratMapping<TIn1, TIn2, TOut>::displayTransformMatrices(const std::string &label) {
		std::cout << label << std::endl;
		std::cout << "Frames exponential SE3 matrices (size: " << m_frames_exponential_se3_vectors.size() << "):" << std::endl;
		for (size_t i = 0; i < m_frames_exponential_se3_vectors.size(); ++i) {
			std::cout << "  Frame[" << i << "]: " << m_frames_exponential_se3_vectors[i] << std::endl;
		}
		std::cout << "Nodes exponential SE3 matrices (size: " << m_nodes_exponential_se3_vectors.size() << "):" << std::endl;
		for (size_t i = 0; i < m_nodes_exponential_se3_vectors.size(); ++i) {
			std::cout << "  Node[" << i << "]: " << m_nodes_exponential_se3_vectors[i] << std::endl;
		}
	}

} // namespace Cosserat::mapping
