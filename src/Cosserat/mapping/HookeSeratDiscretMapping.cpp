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
 ******************************************************************************/
#define SOFA_COSSERAT_CPP_HookeSeratDiscretMapping
#include <Cosserat/mapping/HookeSeratDiscretMapping.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace Cosserat::mapping {
	template class SOFA_COSSERAT_API HookeSeratDiscretMapping<
			sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
	// template class SOFA_COSSERAT_API HookeSeratDiscretMapping<
	// 		sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;

} // namespace Cosserat::mapping

namespace Cosserat {
	// Register in the Factory
	void registerHookeSeratDiscretMapping(sofa::core::ObjectFactory *factory) {
		factory->registerObjects(
				sofa::core::ObjectRegistrationData(
						"This component facilitates the creation of Hooke Serat Discrete Mapping in SOFA simulations. "
						"It takes two mechanical"
						"objects as inputs: the rigid base of the beam (with 6 degrees of freedom) and the local "
						"coordinates of the beam. Using "
						"these inputs, the component computes and outputs the mechanical positions of the beam in "
						"global coordinates. "
						"Like any mapping, it updates the positions and velocities of the outputs based on the inputs. "
						"Additionally, forces applied to the outputs are propagated back to the inputs, ensuring "
						"bidirectional coupling.")
						.add<mapping::HookeSeratDiscretMapping<sofa::defaulttype::Vec3Types,
															   sofa::defaulttype::Rigid3Types,
															   sofa::defaulttype::Rigid3Types>>(true));
		// .add<mapping::HookeSeratDiscretMapping<sofa::defaulttype::Vec6Types,
		// 									   sofa::defaulttype::Rigid3Types,
		// 									   sofa::defaulttype::Rigid3Types>>());
	}
} // namespace Cosserat
