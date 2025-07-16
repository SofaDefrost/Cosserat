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
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <sofa/helper/ColorMap.h>

namespace Cosserat::mapping {

/**
 * @brief Discrete implementation of HookeSeratBaseMapping using liegroups library
 * 
 * This class provides a concrete implementation of the Cosserat rod mapping
 * using the liegroups library for SE(3) operations, with discrete exponential
 * integration along the rod.
 *
 * @tparam TIn1 The first input type for the mapping (strain state)
 * @tparam TIn2 The second input type for the mapping (rigid base)  
 * @tparam TOut The output type for the mapping (frames)
 */
template<class TIn1, class TIn2, class TOut>
class HookeSeratDiscretMapping : public HookeSeratBaseMapping<TIn1, TIn2, TOut> {
public:
    SOFA_CLASS(SOFA_TEMPLATE3(HookeSeratDiscretMapping, TIn1, TIn2, TOut),
               SOFA_TEMPLATE3(HookeSeratBaseMapping, TIn1, TIn2, TOut));

    using In1 = TIn1;
    using In2 = TIn2;
    using Out = TOut;
    using Inherit = HookeSeratBaseMapping<TIn1, TIn2, TOut>;

    // Type aliases from base classes
    using Coord1 = sofa::Coord_t<In1>;
    using Deriv1 = sofa::Deriv_t<In1>;
    using OutCoord = sofa::Coord_t<Out>;
    using OutDeriv = sofa::Deriv_t<Out>;

    //using SectionProperties = typename HookeSeratBaseMapping<TIn1,TIn2,TOut>::SectionProperties;
    //using FrameInfo = typename FrameInfo;
    using SE3Types = sofa::component::cosserat::liegroups::SE3<double>;
    using Vector3 = typename SE3Types::Vector3;
    using Vector6 = typename SE3Types::TangentVector;

public:
    //////////////////////////////////////////////////////////////////////
    /// @name Data Fields
    /// @{
    sofa::Data<int> d_deformationAxis;
    sofa::Data<SReal> d_max;
    sofa::Data<SReal> d_min;
    sofa::Data<SReal> d_radius;
    sofa::Data<bool> d_drawMapBeam;
    sofa::Data<sofa::type::RGBAColor> d_color;
    sofa::Data<sofa::type::vector<int>> d_index;
    sofa::Data<unsigned int> d_baseIndex;
    /// @}
    //////////////////////////////////////////////////////////////////////

public:
    //////////////////////////////////////////////////////////////////////
    /// @name Inherited from BaseObject
    /// @{
    void doBaseCosseratInit() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;
    /// @}
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    /// @name Inherited from Multi2Mapping  
    /// @{
    void apply(const sofa::core::MechanicalParams* mparams,
               const sofa::type::vector<sofa::DataVecCoord_t<Out>*>& dataVecOutPos,
               const sofa::type::vector<const sofa::DataVecCoord_t<In1>*>& dataVecIn1Pos,
               const sofa::type::vector<const sofa::DataVecCoord_t<In2>*>& dataVecIn2Pos) override;

    void applyJ(const sofa::core::MechanicalParams* mparams,
                const sofa::type::vector<sofa::DataVecDeriv_t<Out>*>& dataVecOutVel,
                const sofa::type::vector<const sofa::DataVecDeriv_t<In1>*>& dataVecIn1Vel,
                const sofa::type::vector<const sofa::DataVecDeriv_t<In2>*>& dataVecIn2Vel) override;

    void applyJT(const sofa::core::MechanicalParams* mparams,
                 const sofa::type::vector<sofa::DataVecDeriv_t<In1>*>& dataVecOut1Force,
                 const sofa::type::vector<sofa::DataVecDeriv_t<In2>*>& dataVecOut2RootForce,
                 const sofa::type::vector<const sofa::DataVecDeriv_t<Out>*>& dataVecInForce) override;

    void applyDJT(const sofa::core::MechanicalParams* /*mparams*/, 
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) override {}

    /// Support for constraints
    void applyJT(const sofa::core::ConstraintParams* cparams,
                 const sofa::type::vector<sofa::DataMatrixDeriv_t<In1>*>& dataMatOut1Const,
                 const sofa::type::vector<sofa::DataMatrixDeriv_t<In2>*>& dataMatOut2Const,
                 const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out>*>& dataMatInConst) override;
    /// @}
    //////////////////////////////////////////////////////////////////////

    void computeBBox(const sofa::core::ExecParams* params, bool onlyVisible) override;

protected:
    ////////////////////////// Inherited attributes ////////////////////////////
    /// Bring inherited attributes into the current lookup context
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_sectionProperties;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_frameProperties;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_indicesVectors;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_indicesVectorsDraw;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_beamLengthVectors;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::d_curv_abs_section;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::d_curv_abs_frames;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::d_debug;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_strain_state;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_rigid_base;
    using HookeSeratBaseMapping<TIn1, TIn2, TOut>::m_frames;
    //////////////////////////////////////////////////////////////////////////////

    sofa::helper::ColorMap m_colorMap;

    /**
     * @brief Updates frame transformations using liegroups SE(3) exponential map
     * @param strainState Current strain values
     */
    void updateFrameTransformations(const sofa::type::vector<Coord1>& strainState);

    /**
     * @brief Computes SE(3) exponential at a specific position along a section
     * @param sectionLength Length parameter
     * @param strain Strain vector (3D or 6D depending on TIn1)
     * @return SE(3) transformation
     */
    SE3Types computeSE3Exponential(double sectionLength, const Coord1& strain);

    // Debug display functions
    void displayStrainState(const sofa::type::vector<Coord1>& strainState, const std::string& context = "") const;
    void displayRigidState(const sofa::type::vector<sofa::Coord_t<In2>>& rigidState, const std::string& context = "") const;
    void displayOutputFrames(const sofa::type::vector<OutCoord>& outputFrames, const std::string& context = "") const;
    void displaySectionProperties(const std::string& context = "") const;
    void displayFrameProperties(const std::string& context = "") const;
    void displaySE3Transform(const SE3Types& transform, const std::string& name = "Transform") const;
    void displayMappingState(const std::string& context = "") const;
    void displayVelocities(const sofa::type::vector<Deriv1>& strainVel, 
                          const sofa::type::vector<sofa::Deriv_t<In2>>& baseVel,
                          const sofa::type::vector<OutDeriv>& outputVel,
                          const std::string& context = "") const;

protected:
    HookeSeratDiscretMapping();
    ~HookeSeratDiscretMapping() override = default;
};

#if !defined(SOFA_COSSERAT_CPP_HookeSeratDiscretMapping)
extern template class SOFA_COSSERAT_API HookeSeratDiscretMapping<
    sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
extern template class SOFA_COSSERAT_API HookeSeratDiscretMapping<
    sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
