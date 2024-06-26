/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture                          *
 *                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
 *                                                                             *
 * This library is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this library; if not, write to the Free Software Foundation,     *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
 *******************************************************************************
 *                           Plugin Cosserat    v1.0                           *
 *				                                                              *
 * This plugin is also distributed under the GNU LGPL (Lesser General          *
 * Public License) license with the same conditions than SOFA.                 *
 *                                                                             *
 * Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
 *               Ecole Centrale de Lille)                                      *
 *                                                                             *
 * Contact information: https://project.inria.fr/softrobot/contact/            *
 *                                                                             *
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <sofa/type/vector.h>
#include <sofa/type/Mat.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <Eigen/Dense>
namespace Cosserat
{

    namespace type
    {
        using SE3 =  sofa::type::Matrix4;

        typedef typename sofa::defaulttype::SolidTypes<SReal>::Transform Transform;
        typedef typename sofa::type::vector<SReal> List;

        typedef typename Eigen::Matrix3d RotMat;
        typedef typename Eigen::Matrix<SReal, 6, 1> Vector6d;

        //using TangentTransform = sofa::type::Mat6x6;
        class TangentTransform : public sofa::type::Mat<6, 6, SReal>
        {
        public:
            //TangentTransform(SReal v=0) : sofa::type::Mat6x6(v){}
        };

    }
}

