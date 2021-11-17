//
// Created by younes on 17/11/2021.
//

#include "LegendrePolynomialsMapping.h"

#pragma once
#include <SofaRigid/RigidMapping.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/State.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/helper/io/XspLoader.h>
#include <sofa/helper/io/SphereLoader.h>
#include <sofa/helper/io/Mesh.h>
#include <sofa/helper/decompose.h>
#include <sofa/core/MechanicalParams.h>

namespace sofa::component::mapping {

    template <class TIn, class TOut>
    LegendrePolynomialsMapping<TIn, TOut>::LegendrePolynomialsMapping()
        : Inherit()
        , points(initData(&points, "initialPoints", "Local Coordinates of the points"))
        , index(initData(&index, (unsigned)0, "index", "input DOF index"))

    {
    }

}