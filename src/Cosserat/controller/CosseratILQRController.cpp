/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture                *
 *                 (c) 2006 INRIA, USTL, UJF, CNRS, MGH                       *
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
 ******************************************************************************/
#define SOFA_COSSERAT_CPP_CosseratILQRController

#include <Cosserat/controller/CosseratILQRController.inl>
#include <sofa/core/ObjectFactory.h>

namespace Cosserat::controller {

	using namespace sofa::defaulttype;

	// ── Explicit instantiation ────────────────────────────────────────────────

	template class SOFA_COSSERAT_API CosseratILQRController<
		Vec3Types, Rigid3Types, Rigid3Types>;

	// ── SOFA factory registration ─────────────────────────────────────────────

	void registerCosseratILQRController(sofa::core::ObjectFactory *factory) {
		factory->registerObjects(sofa::core::ObjectRegistrationData(
			"Quasi-static iLQR tip-tracking controller for Cosserat rods.\n"
			"\n"
			"Computes optimal strain corrections at each simulation step to\n"
			"drive the rod tip toward a desired target pose.  Two modes:\n"
			"  mode=0  Gradient descent (fast, may need small stepSize)\n"
			"  mode=1  Gauss-Newton    (quadratic convergence near solution)\n"
			"\n"
			"Requires a link to a Strain2FramesCosseratMapping (mapping=@...).\n"
			"The strain corrections are written back to the strain mechanical\n"
			"state at the start of each AnimateBeginEvent.")
			.add<CosseratILQRController<Vec3Types, Rigid3Types, Rigid3Types>>(true));
	}

} // namespace Cosserat::controller
