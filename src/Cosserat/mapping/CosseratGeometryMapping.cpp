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
#define SOFA_COSSERAT_CPP_CosseratGeometryMapping
#include <Cosserat/mapping/CosseratGeometryMapping.inl>
#include <sofa/defaulttype/VecTypes.h>

namespace Cosserat::mapping {
	// Explicit instantiation of the CosseratBeamGeometry mixin for each TStrain
	// used by the mapping specialisations across the plugin. With Windows DLL
	// exports, the linker needs the mixin's methods materialised in this
	// translation unit. CosseratBeamGeometry<Rigid3> is no longer needed since
	// the option-B refactor removed the old Frames2StrainCosseratMapping that
	// inherited from CosseratGeometryMapping<Rigid3, Rigid3, Vec3/Vec6>.
	template class SOFA_COSSERAT_API CosseratBeamGeometry<sofa::defaulttype::Vec3Types>;
	template class SOFA_COSSERAT_API CosseratBeamGeometry<sofa::defaulttype::Vec6Types>;

	// for Strain2FramesCosseratMapping (Multi2Mapping<Vec3/Vec6, Rigid3, Rigid3>)
	template class SOFA_COSSERAT_API CosseratGeometryMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
															sofa::defaulttype::Rigid3Types>;
	template class SOFA_COSSERAT_API CosseratGeometryMapping<sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types,
															sofa::defaulttype::Rigid3Types>;
} // namespace Cosserat::mapping
