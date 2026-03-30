// This file defines the main pybind11 module for the Cosserat plugin, including bindings for Lie groups and
// PointsManager.

/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture                *
 *                    (c) 2021 INRIA, USTL, UJF, CNRS, MGH                     *
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
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/

#include <pybind11/pybind11.h>
#include "Binding_LieGroups.h"

namespace sofapython3 {

	PYBIND11_MODULE(LieGroups, m) {
		m.doc() = "Cosserat plugin for SOFA, providing Lie group functionalities for Cosserat models.";
		// Add all Lie group bindings
		moduleAddLieGroups(m);
	        // @TODO, add the reste of the bindings here !

	}

} // namespace sofapython3
