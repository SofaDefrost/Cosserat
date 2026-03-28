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
#include <Cosserat/mapping/DifferenceMultiMapping.h>

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/RGBAColor.h>

#include <string>

namespace Cosserat::mapping {
	using sofa::core::objectmodel::BaseContext;
	using sofa::helper::AdvancedTimer;
	using sofa::helper::WriteAccessor;
	using sofa::type::RGBAColor;

	template<class TIn1, class TIn2, class TOut>
	DifferenceMultiMapping<TIn1, TIn2, TOut>::DifferenceMultiMapping() :
		d_direction(initData(&d_direction, "direction", "The list of directions of fix points .\n")),
		d_indices(initData(&d_indices, "indices", "Indices of fixe points of the cable")),
		d_radius(initData(&d_radius, 2.0, "radius", "The size of the cable")),
		d_color(initData(&d_color, sofa::type::Vec4f(1, 0, 0, 1), "color", "The color of the cable")),
		d_drawArrows(initData(&d_drawArrows, false, "drawArrows",
							  "Draw constraint arrows between FEM points and their projections on the cable")),
		d_lastPointIsFixed(initData(&d_lastPointIsFixed, true, "lastPointIsFixed",
									"If true, the last point is treated as a fixed bilateral 3D constraint. "
									"If false (needle mode), all points use proximity-based constraints.")),
		m_fromModel1(nullptr), m_fromModel2(nullptr), m_toModel(nullptr) {}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::initiateTopologies() {
		m_toModel = this->getToModels()[0];
		if (!m_toModel) {
			std::cout << " No output mechanical state found. Consider setting the " << this->toModels.getName()
					  << " attribute." << std::endl;
			return;
		}

		if (!d_direction.isSet())
			msg_warning() << "No direction nor indices is given.";
	}

	// _________________________________________________________________________________________

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::init() {
		Inherit1::init();

		if (this->getFromModels1().empty()) {
			msg_error() << "Error while initializing ; input getFromModels1 not found";
			return;
		}

		if (this->getFromModels2().empty()) {
			msg_error() << "Error while initializing ; output getFromModels2 not found";
			return;
		}

		if (this->getToModels().empty()) {
			msg_error() << "Error while initializing ; output Model not found";
			// return;
		}
		m_fromModel1 = this->getFromModels1()[0];
		m_fromModel2 = this->getFromModels2()[0];
		m_toModel = this->getToModels()[0];

		initiateTopologies();
	}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::bwdInit() {}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::reinit() {}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::reset() {
		reinit();
	}

	// ─────────────────────────────────────────────────────────────────────────────
	// computeProximity / computeNeedleProximity — public wrappers
	// ─────────────────────────────────────────────────────────────────────────────

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::computeProximity(const In1VecCoord &x1, const In2VecCoord &x2) {
		computeProximityImpl(x1, x2, /*fixLastPoint=*/true);
	}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::computeNeedleProximity(const In1VecCoord &x1,
																		  const In2VecCoord &x2) {
		computeProximityImpl(x1, x2, /*fixLastPoint=*/false);
	}

	// ─────────────────────────────────────────────────────────────────────────────
	// computeProximityImpl — common implementation
	//
	// For each FEM point P in x1, finds the closest segment of the cable x2
	// using a true 3D Euclidean distance (closest-point-on-segment), then fills
	// the corresponding Constraint entry.
	//
	// If fixLastPoint==true the last FEM point is handled as a fixed bilateral
	// 3D constraint tied to the last cable node (alpha=1), independently of
	// the proximity search.
	// ─────────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::computeProximityImpl(const In1VecCoord &x1, const In2VecCoord &x2,
																		bool fixLastPoint) {
		m_constraints.clear();

		const size_t szFrom = x1.size();
		const size_t szDst = x2.size();
		const vector<Rigid> &direction = d_direction.getValue();

		for (size_t i = 0; i < szFrom; i++) {
			const Coord2 &P = x1[i];
			Constraint constraint;

			// ── Fixed last-point: bilateral 3D constraint on the last cable node ─
			if (fixLastPoint && i == szFrom - 1) {
				Coord1 dirAxe = x2[szDst - 1] - x2[szDst - 2];
				dirAxe.normalize();

				constraint.eid = static_cast<int>(szDst) - 2;
				constraint.alpha = Real(1.0);
				constraint.proj = x2[szDst - 1];
				constraint.Q = P;
				constraint.dist = (x2[szDst - 1] - P).norm();
				constraint.dirAxe = dirAxe;

				if (!direction.empty()) {
					const Quat _ori = direction[szDst - 1].getOrientation();
					Vec3 _vY = _ori.rotate(Vec3(0., 1., 0.));
					_vY.normalize();
					Vec3 _vZ = _ori.rotate(Vec3(0., 0., 1.));
					_vZ.normalize();
					constraint.t1 = _vY;
					constraint.t2 = _vZ;
				} else {
					// arbitrary normal when no frame is provided
					Deriv1 t2 = cross(constraint.t1, dirAxe);
					t2.normalize();
					constraint.t2 = t2;
				}

				m_constraints.push_back(constraint);
				continue;
			}

			// ── Find segment of x2 closest to P using true 3D distance ──────────
			Real min_dist = std::numeric_limits<Real>::max();
			size_t best_j = 0;

			for (size_t j = 0; j < szDst - 1; j++) {
				const Coord1 Q1 = x2[j];
				const Coord1 Q2 = x2[j + 1];
				const Coord1 seg = Q2 - Q1;
				const Real sq_len = dot(seg, seg);
				if (sq_len < Real(1e-12))
					continue;

				// Clamp parameter t to [0,1] → closest point on segment
				Real t = dot(P - Q1, seg) / sq_len;
				t = std::max(Real(0.0), std::min(Real(1.0), t));

				const Coord1 closest = Q1 + seg * t;
				const Real dist = (P - closest).norm();

				if (dist < min_dist) {
					min_dist = dist;
					best_j = j;
				}
			}

			// ── Compute constraint for the closest segment ────────────────────────
			{
				const Coord1 Q1 = x2[best_j];
				const Coord1 Q2 = x2[best_j + 1];
				Coord1 dirAxe = Q2 - Q1;
				const Real length = dirAxe.norm();

				// Unnormalised projection parameter (0=Q1, 1=Q2)
				const Real fact_v = (length > Real(1e-12)) ? dot(P - Q1, dirAxe) / dot(dirAxe, dirAxe) : Real(0.0);

				constraint.eid = static_cast<int>(best_j);
				constraint.Q1Q2 = length;
				constraint.r2 = fact_v;
				constraint.Q = P;

				dirAxe.normalize();
				constraint.dirAxe = dirAxe;

				// Projection of P onto the (now normalised) axis
				const Real alpha_proj = (P - Q1) * dirAxe;
				const Coord1 proj = Q1 + dirAxe * alpha_proj;
				constraint.proj = proj;
				constraint.dist = (P - proj).norm();

				// alpha: normalised ∈ [0,1], weight of node Q1 (eid) in applyJ/applyJT
				const Real alpha_norm = alpha_proj / length;
				constraint.alpha = (alpha_norm < Real(1e-8)) ? Real(1.0) : (Real(1.0) - alpha_norm);

				// Frame tangent direction from the current segment's rigid frame
				if (!direction.empty() && best_j < direction.size()) {
					const Quat ori = direction[best_j].getOrientation();
					Vec3 vY = Vec3(0., 1., 0.);
					vY = ori.rotate(vY);
					vY.normalize();
					constraint.t1 = vY;
				} else {
					// Fallback: use violation vector, or arbitrary if P is on the axis
					Deriv1 t1 = P - proj;
					if (t1.norm() < Real(1e-8))
						t1 = Deriv1(0., 1., 0.);
					t1.normalize();
					constraint.t1 = t1;
				}

				Deriv1 t2 = cross(constraint.t1, dirAxe);
				t2.normalize();
				constraint.t2 = t2;
			}

			m_constraints.push_back(constraint);
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams * /* mparams */,
														 const vector<OutDataVecCoord *> &dataVecOutPos,
														 const vector<const In1DataVecCoord *> &dataVecIn1Pos,
														 const vector<const In2DataVecCoord *> &dataVecIn2Pos) {

		if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		const In1VecCoord &in1 = dataVecIn1Pos[0]->getValue();
		const In2VecCoord &in2 = dataVecIn2Pos[0]->getValue();

		OutVecCoord &out = *dataVecOutPos[0]->beginEdit();

		if (d_lastPointIsFixed.getValue()) {
			computeProximity(in1, in2);

			size_t sz = m_constraints.size();
			out.resize(sz);

			for (unsigned int i = 0; i < sz; i++) {
				Constraint &c = m_constraints[i];
				if (i < sz - 1) {
					out[i][0] = 0.0;
					out[i][1] = c.t1 * (in1[i] - c.proj); // c.dist;
					out[i][2] = c.t2 * (in1[i] - c.proj); // 0.0
				} else {
					out[sz - 1][0] = c.dirAxe * (in1[in1.size() - 1] - in2[in2.size() - 1]);
					out[sz - 1][1] = c.t1 * (in1[in1.size() - 1] - in2[in2.size() - 1]);
					out[sz - 1][2] = c.t2 * (in1[in1.size() - 1] - in2[in2.size() - 1]);
				}
			}
		} else {
			computeNeedleProximity(in1, in2);

			size_t sz = m_constraints.size();
			out.resize(sz);

			for (unsigned int i = 0; i < sz; i++) {
				Constraint &c = m_constraints[i];
				out[i][0] = 0.0;
				out[i][1] = c.t1 * (in1[i] - c.proj); // c.dist;
				out[i][2] = c.t2 * (in1[i] - c.proj); // 0.0
			}
		}
		dataVecOutPos[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams * /* mparams */,
														  const vector<OutDataVecDeriv *> &dataVecOutVel,
														  const vector<const In1DataVecDeriv *> &dataVecIn1Vel,
														  const vector<const In2DataVecDeriv *> &dataVecIn2Vel) {
		if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;
		const In1VecDeriv &in1 = dataVecIn1Vel[0]->getValue();
		const In2VecDeriv &in2 = dataVecIn2Vel[0]->getValue();
		OutVecDeriv &outVel = *dataVecOutVel[0]->beginEdit();

		size_t sz = m_constraints.size();
		outVel.resize(sz);

		if (d_lastPointIsFixed.getValue()) {
			for (size_t i = 0; i < sz; i++) {
				Constraint &c = m_constraints[i];
				int ei1 = c.eid;
				int ei2 = c.eid + 1;
				if (i < sz - 1) {
					Real v0 = c.dirAxe * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
					Real v1 = c.t1 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
					Real v2 = c.t2 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
					outVel[i] = OutDeriv(v0, v1, v2);
				} else {
					Real v0 = c.dirAxe * (in1[i] - in2[ei2]);
					Real v1 = c.t1 * (in1[i] - in2[ei2]);
					Real v2 = c.t2 * (in1[i] - in2[ei2]);
					outVel[i] = OutDeriv(v0, v1, v2);
				}
			}
		} else {
			for (size_t i = 0; i < sz; i++) {
				Constraint &c = m_constraints[i];

				int ei1 = c.eid;
				int ei2 = c.eid + 1;
				Real v0 = c.dirAxe * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
				Real v1 = c.t1 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
				Real v2 = c.t2 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);

				outVel[i] = OutDeriv(v0, v1, v2);
			}
		}
		dataVecOutVel[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::MechanicalParams * /*mparams*/,
														   const vector<In1DataVecDeriv *> &dataVecOut1Force,
														   const vector<In2DataVecDeriv *> &dataVecOut2Force,
														   const vector<const OutDataVecDeriv *> &dataVecInForce) {
		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		const OutVecDeriv &in = dataVecInForce[0]->getValue();
		if (in.empty())
			return;

		In1VecDeriv &out1 = *dataVecOut1Force[0]->beginEdit();
		In2VecDeriv &out2 = *dataVecOut2Force[0]->beginEdit();
		// Compute output forces
		size_t sz = m_constraints.size();
		if (d_lastPointIsFixed.getValue()) {
			for (size_t i = 0; i < sz; i++) {
				Constraint &c = m_constraints[i];
				int ei1 = c.eid;
				int ei2 = c.eid + 1;
				OutDeriv f = in[i];
				if (i < sz - 1) {
					Deriv2 f1 = (f[0] * c.dirAxe) + (f[1] * c.t1) + (f[2] * c.t2);
					Deriv1 f2_1 = (c.alpha * f[0] * c.dirAxe) + (c.alpha * f[1] * c.t1) + (c.alpha * f[2] * c.t2);
					Deriv1 f2_2 = ((1 - c.alpha) * f[0] * c.dirAxe) + ((1 - c.alpha) * f[1] * c.t1) +
								  ((1 - c.alpha) * f[2] * c.t2);
					out1[i] += f1;
					out2[ei1] -= f2_1;
					out2[ei2] -= f2_2;
				} else {
					Deriv2 f1 = (f[0] * c.dirAxe) + (f[1] * c.t1) + (f[2] * c.t2);
					out1[i] += f1;
					out2[ei2] -= f1;
				}
			}
		} else {
			for (size_t i = 0; i < sz; i++) {
				Constraint &c = m_constraints[i];
				int ei1 = c.eid;
				int ei2 = c.eid + 1;
				OutDeriv f = in[i];
				Deriv2 f1 = (f[0] * c.dirAxe) + (f[1] * c.t1) + (f[2] * c.t2);
				Deriv1 f2_1 = (c.alpha * f[0] * c.dirAxe) + (c.alpha * f[1] * c.t1) + (c.alpha * f[2] * c.t2);
				Deriv1 f2_2 = ((1 - c.alpha) * f[0] * c.dirAxe) + ((1 - c.alpha) * f[1] * c.t1) +
							  ((1 - c.alpha) * f[2] * c.t2);
				out1[i] += f1;
				out2[ei1] -= f2_1;
				out2[ei2] -= f2_2;
			}
		}
		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();
	}

	//___________________________________________________________________________
	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::ConstraintParams * /*cparams*/,
														   const vector<In1DataMatrixDeriv *> &dataMatOut1Const,
														   const vector<In2DataMatrixDeriv *> &dataMatOut2Const,
														   const vector<const OutDataMatrixDeriv *> &dataMatInConst) {
		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;

		// We need only one input In model and input Root model (if present)
		In1MatrixDeriv &out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the FEM cable points
		In2MatrixDeriv &out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the frames cable points
		const OutMatrixDeriv &in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped point
		const In1DataVecCoord *x1fromData = m_fromModel1->read(sofa::core::vec_id::read_access::position);
		const In1VecCoord x1from = x1fromData->getValue();

		typename OutMatrixDeriv::RowConstIterator rowIt = in.begin();
		typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

		for (rowIt = in.begin(); rowIt != rowItEnd; ++rowIt) {
			typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
			typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

			// Creates a constraints if the input constraint is not empty.
			if (colIt == colItEnd) {
				continue;
			}
			typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
			typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());

			if (d_lastPointIsFixed.getValue()) {
				if ((rowIt.index() / 2) < (int) (x1from.size() - 1)) {
					while (colIt != colItEnd) {
						int childIndex = colIt.index();
						Constraint c = m_constraints[childIndex];
						const OutDeriv h = colIt.val();
						int indexBeam = c.eid;

						Deriv2 h1 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);
						Deriv1 h2_1 = (c.alpha * h[0] * c.dirAxe) + (c.alpha * h[1] * c.t1) + (c.alpha * h[2] * c.t2);
						Deriv1 h2_2 = ((1.0 - c.alpha) * h[0] * c.dirAxe) + ((1.0 - c.alpha) * h[1] * c.t1) +
									  ((1.0 - c.alpha) * h[2] * c.t2);

						o1.addCol(childIndex, h1);
						o2.addCol(indexBeam, -h2_1);
						o2.addCol(indexBeam + 1, -h2_2);

						colIt++;
					}
				} else {
					while (colIt != colItEnd) {
						int childIndex = colIt.index();
						Constraint c = m_constraints[childIndex];
						const OutDeriv h = colIt.val();
						int indexBeam = c.eid;

						Deriv2 h1 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);
						Deriv1 h2 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);

						o1.addCol(childIndex, h1);
						o2.addCol(indexBeam + 1, -h2);
						colIt++;
					}
				}
			} else {
				while (colIt != colItEnd) {
					int childIndex = colIt.index();
					Constraint c = m_constraints[childIndex];
					const OutDeriv h = colIt.val();
					int indexBeam = c.eid;

					Deriv2 h1 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);
					Deriv1 h2_1 = (c.alpha * h[0] * c.dirAxe) + (c.alpha * h[1] * c.t1) + (c.alpha * h[2] * c.t2);
					Deriv1 h2_2 = ((1.0 - c.alpha) * h[0] * c.dirAxe) + ((1.0 - c.alpha) * h[1] * c.t1) +
								  ((1.0 - c.alpha) * h[2] * c.t2);

					o1.addCol(childIndex, h1);
					o2.addCol(indexBeam, -h2_1);
					o2.addCol(indexBeam + 1, -h2_2);

					colIt++;
				}
			}
		}
		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void DifferenceMultiMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		/// draw cable
		if (!vparams->displayFlags().getShowInteractionForceFields())
			return;

		typedef sofa::type::RGBAColor RGBAColor;
		vparams->drawTool()->saveLastState();
		vparams->drawTool()->disableLighting();

		std::vector<Vec3> vertices;
		RGBAColor color = RGBAColor::magenta();

		if (d_drawArrows.getValue() && d_lastPointIsFixed.getValue()) {
			for (size_t i = 0; i < m_constraints.size(); i++) {
				color = RGBAColor::green();
				vertices.push_back(m_constraints[i].proj);
				vertices.push_back(m_constraints[i].Q);
				vparams->drawTool()->drawLines(vertices, 4.0, color);
				if (i == (m_constraints.size() - 1)) {
					Coord2 P1 = m_constraints[i].Q;
					Real radius_arrow = 0.30;
					Coord2 x = m_constraints[i].dirAxe * 5.0;
					Coord2 y = m_constraints[i].t1 * 5.0;
					Coord2 z = m_constraints[i].t2 * 5.0;

					vparams->drawTool()->drawArrow(P1, P1 + x, radius_arrow, RGBAColor::red());
					vparams->drawTool()->drawArrow(P1, P1 + y, radius_arrow, RGBAColor::green());
					vparams->drawTool()->drawArrow(P1, P1 + z, radius_arrow, RGBAColor::blue());
				} else {
					Coord2 P1 = m_constraints[i].Q;
					Real radius_arrow = 0.30;
					Coord2 y = m_constraints[i].t1 * 5.0;
					Coord2 z = m_constraints[i].t2 * 5.0;
					vparams->drawTool()->drawArrow(P1, P1 + y, radius_arrow, RGBAColor::blue());
					vparams->drawTool()->drawArrow(P1, P1 + z, radius_arrow, RGBAColor::blue());
				}
			}
			const In1DataVecDeriv *xDestData = m_fromModel1->read(sofa::core::vec_id::read_access::position);
			const In1VecCoord &fromPos = xDestData[0].getValue();
			vparams->drawTool()->draw3DText_Indices(fromPos, 6, RGBAColor(0.0, 1.0, 0.0, 1.0));
		}
		vparams->drawTool()->restoreLastState();
	}

} // namespace Cosserat::mapping
