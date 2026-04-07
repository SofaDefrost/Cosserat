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
#include <Cosserat/mapping/DiscreteDynamicCosseratMapping.h>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/type/Quat.h>


namespace Cosserat::mapping {
	using sofa::core::objectmodel::BaseContext;
	using sofa::helper::AdvancedTimer;
	using sofa::helper::WriteAccessor;


	template<class TIn1, class TIn2, class TOut>
	DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::DiscreteDynamicCosseratMapping() = default;

	template<class TIn1, class TIn2, class TOut>
	void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::doBaseCosseratInit() {
		for (size_t i = 0; i < 3; i++)
			m_matrixBi[i][i] = 1.0;
	}


	template<class TIn1, class TIn2, class TOut>
	void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams * /* mparams */,
																 const vector<OutDataVecCoord *> &dataVecOutPos,
																 const vector<const In1DataVecCoord *> &dataVecIn1Pos,
																 const vector<const In2DataVecCoord *> &dataVecIn2Pos) {

		if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		/// Do Apply
		// We need only one input In model and input Root model (if present)
		const In1VecCoord &in1 = dataVecIn1Pos[0]->getValue();
		const In2VecCoord &in2 = dataVecIn2Pos[0]->getValue();

		size_t sz = d_curv_abs_frames.getValue().size();
		OutVecCoord out = sofa::helper::getWriteAccessor(*dataVecOutPos[0]);
		out.resize(sz);

		// update the Exponential Matrices according to new deformation
		this->updateExponentialSE3(in1);

		Transform frame0 = Transform(In2::getCPos(in2[0]), In2::getCRot(in2[0]));
		for (unsigned int i = 0; i < sz; i++) {
			Transform frame = frame0;
			for (unsigned int u = 0; u < m_indicesVectors[i]; u++) {
				frame *= m_nodesExponentialSE3Vectors[u];
			}
			frame *= m_framesExponentialSE3Vectors[i];

			Vec3 v = frame.getOrigin();
			sofa::type::Quat q = frame.getOrientation();
			out[i] = OutCoord(v, q);
		}

		m_indexInput = 0;
	}


	template<class TIn1, class TIn2, class TOut>
	void
	DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams * /* mparams */,
															 const vector<OutDataVecDeriv *> &dataVecOutVel,
															 const vector<const In1DataVecDeriv *> &dataVecIn1Vel,
															 const vector<const In2DataVecDeriv *> &dataVecIn2Vel) {

		if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;
		const In1VecDeriv &in1 = dataVecIn1Vel[0]->getValue();
		const In2VecDeriv &in2_vecDeriv = dataVecIn2Vel[0]->getValue();
		OutVecDeriv &outVel = *dataVecOutVel[0]->beginEdit();

		sofa::helper::ReadAccessor<Data<vector<double>>> curv_abs_input =
				d_curv_abs_section; // This is the vector of X in the paper
		sofa::helper::ReadAccessor<Data<vector<double>>> curv_abs_output = d_curv_abs_frames;
		sofa::helper::ReadAccessor<Data<bool>> debug = d_debug;
		m_frameJacobienVector.clear();
		m_frameJacobienDotVector.clear();

		// Compute the tangent Exponential SE3 vectors
		const In1VecCoord &inDeform = m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();
		this->updateTangExpSE3(inDeform);

		// Get base velocity as input this is also called eta
		m_nodesVelocityVectors.clear();
		Deriv2 _baseVelocity;
		if (!in2_vecDeriv.empty())
			_baseVelocity = in2_vecDeriv[0];

		// convert to Vec6
		Vec6 baseVelocity;
		for (size_t u = 0; u < 6; u++) {
			baseVelocity[u] = _baseVelocity[u];
		}

		// Apply the local transform i.e from SOFA frame to Frederico frame
		const In2VecCoord &xfrom2Data = m_rigid_base->read(sofa::core::vec_id::read_access::position)->getValue();
		Transform Tinverse = Transform(xfrom2Data[0].getCenter(), xfrom2Data[0].getOrientation()).inversed();
		Mat6x6 P = this->buildProjector(Tinverse);
		m_nodeAdjointVectors.clear();

		Vec6 baseLocalVelocity = P * baseVelocity;
		m_nodesVelocityVectors.push_back(baseLocalVelocity);
		if (debug)
			std::cout << "Base local Velocity :" << baseLocalVelocity << std::endl;

		// Compute velocity at nodes
		for (size_t i = 1; i < curv_abs_input.size(); i++) {
			Transform t = m_nodesExponentialSE3Vectors[i].inversed();
			TangentTransform Adjoint;
			Adjoint.clear();
			this->computeAdjoint(t, Adjoint);
			// Add this line because need for the jacobian computation
			m_nodeAdjointVectors.push_back(Adjoint);

			// Compute velocity (eta) at node i != 0 eq.(13) paper
			Vec6 Xi_dot = Vec6(in1[i - 1], Vec3(0.0, 0.0, 0.0));
			Vec6 temp = Adjoint * (m_nodesVelocityVectors[i - 1] + m_nodesTangExpVectors[i] * Xi_dot);
			m_nodesVelocityVectors.push_back(temp);
			if (debug)
				std::cout << "Node velocity : " << i << " = " << temp << std::endl;
		}

		const OutVecCoord &out = m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		size_t sz = curv_abs_output.size();
		outVel.resize(sz);
		for (size_t i = 0; i < sz; i++) {
			Transform transform = m_framesExponentialSE3Vectors[i].inversed();
			TangentTransform tangentTransform;
			tangentTransform.clear();
			this->computeAdjoint(transform, tangentTransform);

			Vec6 Xi_dot = Vec6(in1[m_indicesVectors[i] - 1], Vec3(0.0, 0.0, 0.0));
			// eta
			Vec6 etaFrame = tangentTransform *
							(m_nodesVelocityVectors[m_indicesVectors[i] - 1] + m_framesTangExpVectors[i] * Xi_dot);

			// Compute here jacobien and jacobien_dot
			vector<Mat6x3> J_i, J_dot_i;
			computeJ_Jdot_i(tangentTransform, i, J_i, etaFrame, J_dot_i);
			m_frameJacobienVector.push_back(J_i);
			m_frameJacobienVector.push_back(J_dot_i);

			// Convert from Federico node to Sofa node
			Transform _T = Transform(out[i].getCenter(), out[i].getOrientation());
			Mat6x6 _P = this->buildProjector(_T);
			// std::cout<< "Eta local : "<< eta << std::endl;

			outVel[i] = _P * etaFrame;
			if (debug)
				std::cout << "Frame velocity : " << i << " = " << etaFrame << std::endl;
		}
		//    std::cout << "Inside the apply J, outVel after computation  :  "<< outVel << std::endl;
		dataVecOutVel[0]->endEdit();
		m_indexInput = 0;
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::computeJ_Jdot_i(const Mat6x6 &Adjoint, size_t frameId,
																		   vector<Mat6x3> &J_i, const Vec6 &etaFrame,
																		   vector<Mat6x3> &J_dot_i) {
		size_t sz = d_curv_abs_section.getValue().size();
		Mat6x3 Si;
		Mat6x6 M;

		Mat6x3 Si_dot;
		Mat6x6 adj_eta; // to be computed

		bool reachNode = false;
		for (unsigned int i = 1; i < sz; i++) {
			M = Adjoint;
			unsigned int u = m_indicesVectors[frameId];
			if (i < u) {
				for (int j = u; j > 0; j--) {
					M = M * m_nodeAdjointVectors[j];
				}
				Mat6x6 temp = M * m_nodesTangExpVectors[u];
				Si = temp * m_matrixBi;
				J_i.push_back(Si);

				Vec6 etaNode = m_nodesVelocityVectors[i];
				this->computeAdjoint(etaNode, adj_eta);
				Si_dot = temp * adj_eta * m_matrixBi;
				J_dot_i.push_back(Si_dot);
			} else {
				if (!reachNode) {
					Mat6x6 temp = M * m_framesTangExpVectors[frameId];
					Si = temp * m_matrixBi;
					J_i.push_back(Si);

					this->computeAdjoint(etaFrame, adj_eta);
					Si_dot = temp * adj_eta * m_matrixBi;
					J_dot_i.push_back(Si_dot);
					reachNode = true;
				} else {
					Si.clear();
					J_i.push_back(Si);
					Si_dot.clear();
					J_dot_i.push_back(Si);
				}
			}
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/, const vector<In1DataVecDeriv *> &dataVecOut1Force,
			const vector<In2DataVecDeriv *> &dataVecOut2Force, const vector<const OutDataVecDeriv *> &dataVecInForce) {

		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		const OutVecDeriv &in = dataVecInForce[0]->getValue();

		In1VecDeriv &out1 = *dataVecOut1Force[0]->beginEdit();
		In2VecDeriv &out2 = *dataVecOut2Force[0]->beginEdit();

		// Maybe need, in case the apply funcion is not call this must be call before
		const OutVecCoord &frame = m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		const In1DataVecCoord *x1fromData = m_strain_state->read(sofa::core::vec_id::read_access::position);
		const In1VecCoord x1from = x1fromData->getValue();
		vector<Vec6> cartesianForcesLocalFrame;
		cartesianForcesLocalFrame.clear();

		out1.resize(x1from.size());

		// Step 1: Convert Cartesian environmental forces (input) from SOFA's global space to the beam's local reference
		// frame
		for (size_t var = 0; var < in.size(); ++var) {
			Vec6 forceGlobal;
			for (unsigned j = 0; j < 6; j++)
				forceGlobal[j] = in[var][j];

			Transform frameTransform = Transform(frame[var].getCenter(), frame[var].getOrientation());
			Mat6x6 projectorTransposed = (this->buildProjector(frameTransform));
			projectorTransposed.transpose();

			Vec6 localForce = projectorTransposed * forceGlobal;
			cartesianForcesLocalFrame.push_back(localForce);
		}

		// Step 2: Compute output forces projected into the strain parameter space
		size_t sz = m_indicesVectors.size();

		unsigned int index = m_indicesVectors[sz - 1];
		m_totalBeamForceVectors.clear();
		m_totalBeamForceVectors.resize(sz);

		Vec6 F_tot;
		F_tot.clear();
		m_totalBeamForceVectors.push_back(F_tot);

		// matrix Bi (transposed) extracts the 3 relevant strain components (e.g. elongation/bending) from the full 6D
		// Cartesian twist/wrench
		Mat3x6 strainSelectorBi;
		strainSelectorBi.clear();
		for (unsigned int k = 0; k < 3; k++)
			strainSelectorBi[k][k] = 1.0;

		// Backward propagation of forces from the tip to the root
		for (unsigned int s = sz; s--;) {
			Mat6x6 adjointFrame;
			this->computeCoAdjoint(m_framesExponentialSE3Vectors[s], adjointFrame);

			// Transfer Cartesian force from the beam's geometric frame to its underlying computational node
			Vec6 nodeCartesianForce = adjointFrame * cartesianForcesLocalFrame[s];

			Mat6x6 transposedTangent = m_framesTangExpVectors[s];
			transposedTangent.transpose();

			// Compute force in the strain subspace (reduced coordinates)
			Vec3 strainForceAtFrame = strainSelectorBi * transposedTangent * nodeCartesianForce;

			// If segments bypass structural nodes, correctly propagate the accumulated force downward
			while (index != m_indicesVectors[s]) {
				index--;
				Mat6x6 adjointNode;
				this->computeCoAdjoint(m_nodesExponentialSE3Vectors[index], adjointNode);

				// Translate accumulated force 'F_tot' to the adjacent inferior node
				F_tot = adjointNode * F_tot;

				Mat6x6 transposedNodeTangent = m_nodesTangExpVectors[index];
				transposedNodeTangent.transpose();

				// Map accumulated Cartesian forces to the strain space
				Vec3 strainForceAtNode = strainSelectorBi * transposedNodeTangent * F_tot;
				out1[index - 1] += strainForceAtNode;
			}
			if (d_debug.getValue())
				std::cout << "f at s =" << s << " and index" << index << " is : " << strainForceAtFrame << std::endl;

			// Accumulate nodal Cartesian forces
			F_tot += nodeCartesianForce;
			out1[m_indicesVectors[s] - 1] += strainForceAtFrame;
		}

		// Finally, propagate the total resultant force to the Rigid base of the entire Cosserat structure
		Transform frame0 = Transform(frame[0].getCenter(), frame[0].getOrientation());
		Mat6x6 baseProjectorM = this->buildProjector(frame0);
		out2[0] += baseProjectorM * F_tot;

		if (d_debug.getValue()) {
			std::cout << "Node forces " << out1 << std::endl;
			std::cout << "base Force: " << out2[0] << std::endl;
		}

		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();
	}

	//___________________________________________________________________________
	template<class TIn1, class TIn2, class TOut>
	void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::ConstraintParams * /*cparams*/, const vector<In1DataMatrixDeriv *> &dataMatOut1Const,
			const vector<In2DataMatrixDeriv *> &dataMatOut2Const,
			const vector<const OutDataMatrixDeriv *> &dataMatInConst) {
		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;

		// We need only one input In model and input Root model (if present)
		In1MatrixDeriv &out1 =
				*dataMatOut1Const[0]->beginEdit(); // constraints on the strain space (reduced coordinate)
		In2MatrixDeriv &out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the reference frame (base frame)
		const OutMatrixDeriv &in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

		const OutVecCoord &frame = m_global_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		const In1DataVecCoord *x1fromData = m_strain_state->read(sofa::core::vec_id::read_access::position);
		const In1VecCoord x1from = x1fromData->getValue();
		sofa::helper::ReadAccessor<Data<bool>> debug = d_debug;

		// matrix Bi (transposed) extracts the relevant DOFs from the 6D Cartesian twist mapping
		Mat3x6 strainSelectorBi;
		strainSelectorBi.clear();
		for (unsigned int k = 0; k < 3; k++)
			strainSelectorBi[k][k] = 1.0;

		typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

		for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt) {
			if (debug) {
				std::cout << "************* Apply JT (MatrixDeriv) iteration on line ";
				std::cout << rowIt.index();
				std::cout << "*************  " << std::endl;
			}
			typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
			typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

			// Validates constraints existence
			if (colIt == colItEnd) {
				if (debug) {
					std::cout << "no column for this constraint" << std::endl;
				}
				continue;
			}
			typename In1MatrixDeriv::RowIterator o1 =
					out1.writeLine(rowIt.index()); // Store constraint projection in strain space
			typename In2MatrixDeriv::RowIterator o2 =
					out2.writeLine(rowIt.index()); // Store constraint projection at root

			size_t maxBeamIndex = 0;
			if (!m_indicesVectors.empty())
				maxBeamIndex = m_indicesVectors.back();

			// Accumulator array to securely aggregate overlapping contact forces per-node in O(N) instead of sorting
			// tuples
			std::vector<Vec6> accumulatedCartesianForces(maxBeamIndex + 1);

			while (colIt != colItEnd) {
				int childIndex = colIt.index();

				const OutDeriv valueConst_ = colIt.val();
				Vec6 constraintDirectionGlobal;
				for (unsigned j = 0; j < 6; j++)
					constraintDirectionGlobal[j] = valueConst_[j];

				int indexBeam = m_indicesVectors[childIndex];

				// Step 1: Transfer environmental constraints to the beam's geometric reference
				Transform frameTransform = Transform(frame[childIndex].getCenter(), frame[childIndex].getOrientation());
				Mat6x6 projectorTransposed = (this->buildProjector(frameTransform));
				projectorTransposed.transpose();

				Mat6x6 adjointFrame;
				this->computeCoAdjoint(m_framesExponentialSE3Vectors[childIndex], adjointFrame);

				Mat6x6 transposedTangentFrame = m_framesTangExpVectors[childIndex];
				transposedTangentFrame.transpose();

				// Evaluate Cartesian force dynamically localized to the beam
				Vec6 localConstraintForce = adjointFrame * projectorTransposed * constraintDirectionGlobal;

				// Evaluate corresponding deformation impulse in the strain space
				Vec3 strainConstraintForce = strainSelectorBi * transposedTangentFrame * localConstraintForce;

				if (indexBeam >= 1)
					o1.addCol(indexBeam - 1, strainConstraintForce);

				// Temporarily buffer forces at nodes for rapid backward traversal
				if (indexBeam <= maxBeamIndex)
					accumulatedCartesianForces[indexBeam] += localConstraintForce;

				colIt++;
			}

			// Step 2: Propagate dynamics down to the root structure securely in O(N) linear time (backward accumulation
			// mechanism)
			Vec6 totalBackwardForce;
			totalBackwardForce.clear();
			for (int i = maxBeamIndex; i > 0; i--) {
				totalBackwardForce += accumulatedCartesianForces[i];

				Mat6x6 adjointNode;
				this->computeCoAdjoint(m_nodesExponentialSE3Vectors[i - 1], adjointNode);

				// Transfer the total impulse force downwards continuously
				totalBackwardForce = adjointNode * totalBackwardForce;

				Mat6x6 transposedTangentNode = m_nodesTangExpVectors[i - 1];
				transposedTangentNode.transpose();

				// Re-map the cascading result to affect actual physical strain variables (reduce degrees)
				Vec3 strainForceAtNode = strainSelectorBi * transposedTangentNode * totalBackwardForce;

				if (i > 1) {
					o1.addCol(i - 2, strainForceAtNode);
				}
			}

			totalBackwardForce += accumulatedCartesianForces[0];

			// Finalize transferring constraints onto the global Rigid body (base)
			Transform frame0 = Transform(frame[0].getCenter(), frame[0].getOrientation());
			Mat6x6 baseProjectorM = this->buildProjector(frame0);

			Vec6 constraintForceAtBase = baseProjectorM * totalBackwardForce;
			o2.addCol(0, constraintForceAtBase);
		}

		////// END ARTICULATION SYSTEM MAPPING
		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		if (!vparams->displayFlags().getShowMappings())
			return;
	}

} // namespace Cosserat::mapping
