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
*				                                              *
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
#include <Cosserat/forcefield/BaseBeamHookeLawForceField.h>

#include <sofa/linearalgebra/FullVector.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/OptionsGroup.h>

// Standard includes
#include <iostream>
#include <algorithm>
#include <ctime>

// Using declarations for common types
using sofa::core::behavior::MechanicalState;
using sofa::core::objectmodel::BaseContext;
using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::core::VecCoordId;
using std::cout;
using std::endl;

namespace sofa::component::forcefield
{

using sofa::core::behavior::MultiMatrixAccessor;
using sofa::core::behavior::BaseMechanicalState;
using sofa::helper::WriteAccessor;
/**
 * @brief Get the current frame transformation for a specific cross-section
 * 
 * @param index Cross-section index
 * @return Rigid transformation (SE3) representing the cross-section frame
 */
template<typename DataTypes>
typename BaseBeamHookeLawForceField<DataTypes>::SE3Type 
BaseBeamHookeLawForceField<DataTypes>::getFrame(size_t index) const
{
    assert(index < m_currentFrames.size());
    return m_currentFrames[index];
}

/**
 * @brief Get the reference frame transformation for a specific cross-section
 * 
 * @param index Cross-section index
 * @return Rigid transformation (SE3) representing the reference cross-section frame
 */
template<typename DataTypes>
typename BaseBeamHookeLawForceField<DataTypes>::SE3Type 
BaseBeamHookeLawForceField<DataTypes>::getReferenceFrame(size_t index) const
{
    assert(index < m_referenceFrames.size());
    return m_referenceFrames[index];
}

/**
 * @brief Get the relative rotation between consecutive cross-sections
 * 
 * @param index Segment index (between cross-sections index and index+1)
 * @return Rotation (SO3) representing the relative orientation
 */
template<typename DataTypes>
typename BaseBeamHookeLawForceField<DataTypes>::SO3Type 
BaseBeamHookeLawForceField<DataTypes>::getRelativeRotation(size_t index) const
{
    assert(index < m_relativeRotations.size());
    return m_relativeRotations[index];
}

} // namespace sofa::component::forcefield
 * (different properties for different segments).
 */

/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS                                                 *
******************************************************************************/
ttemplate<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::init()
{
    Inherit1::init();
    if (!validateMaterialParameters()) {
        msg_error("BaseBeamHookeLawForceField") << "Invalid material parameters detected during initialization.";
        return;
    }
    
    if (!validateSectionParameters()) {
        msg_error("BaseBeamHookeLawForceField") << "Invalid cross-section parameters detected during initialization.";
        return;
    }
    
    computeSectionProperties();
    
    // Initialize reference frames
    const VecCoord& x0 = this->mstate->read(sofa::core::vec_id::read_access::restPosition)->getValue();
    const size_t numNodes = x0.size();
    
    m_referenceFrames.resize(numNodes);
    m_currentFrames.resize(numNodes);
    m_relativeRotations.resize(numNodes > 0 ? numNodes - 1 : 0);
    
    // Set up initial reference frames
    for (size_t i = 0; i < numNodes; ++i) {
        // Extract position and rotation from x0[i] 
        Vector3 position = getPosition(x0[i]);
        SO3Type rotation = getRotation(x0[i]);
        
        // Create an SE3 transformation
        m_referenceFrames[i] = SE3Type(rotation, position);
        m_currentFrames[i] = m_referenceFrames[i]; // Initially, current frames match reference frames
    }
    
    // Compute initial relative rotations
    computeRelativeRotations();
    
    m_initialized = true;
}
// For variant section
                localForce = -(m_K_sectionList[i] * strains[i]) * lengths[i];
            }
            
            // Transform to global frame and apply to nodes
            Deriv globalForce = transformWrenchToGlobal(localForce, i);
            
            // Apply forces to both nodes of the segment
            forces[i] += globalForce;
            forces[i+1] -= globalForce;
        }
        
        d_f.endEdit();
    }
    else {
        // Call the original specialized force calculation method
        if (!d_variantSections.getValue()) {
            addForceUniformSection(d_f, d_x, d_x0, lengths);
        } else {
            addForceVariantSection(d_f, d_x, d_x0, lengths);
        }
    }
}
* INITIALIZATION METHODS                                                     *
******************************************************************************/

/**
 * @brief Initialize the force field, setting up parameters and cross-section properties
 * 
 * This method is called during initialization to set up the force field.
 * It validates parameters and calculates cross-section properties.
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::init()
{
    Inherit1::init();
    if (!validateMaterialParameters()) {
        msg_error("BaseBeamHookeLawForceField") << "Invalid material parameters detected during initialization.";
        return;
    }
    
    if (!validateSectionParameters()) {
        msg_error("BaseBeamHookeLawForceField") << "Invalid cross-section parameters detected during initialization.";
        return;
    }
    
    computeSectionProperties();
    m_initialized = true;
}
/**
 * @brief Main method for force computation
 * 
 * This method dispatches the force computation to either addForceUniformSection or
 * addForceVariantSection based on the beam configuration.
 * 
 * @param mparams Mechanical parameters
 * @param d_f Force vector
 * @param d_x Position vector
 * @param d_v Velocity vector (unused)
 */

/**
 * @brief Validates input data for force computation
 * 
 * Checks that the input data is valid for force computation, including:
 * - The mechanical state exists
 * - Position vectors are not empty and have matching sizes
 * - The length vector has the correct size
 * - For variant sections, the section lists have the correct sizes
 * 
 * @param d_f Force vector
 * @param d_x Position vector
 * @param d_x0 Rest position vector
 * @return true if data is valid, false otherwise
 */
template<typename DataTypes>
bool BaseBeamHookeLawForceField<DataTypes>::validateInputData(const DataVecDeriv& d_f, 
                                                             const DataVecCoord& d_x, 
                                                             const DataVecCoord& d_x0) const
{
    if (!this->getMState()) {
        msg_info("BaseBeamHookeLawForceField") << "No Mechanical State found, no force will be computed..." << "\n";
        return false;
    }
    
    const VecCoord& x = d_x.getValue();
    const VecCoord& x0 = d_x0.getValue();
    
    if (x.empty() || x0.empty()) {
        msg_warning("BaseBeamHookeLawForceField") << "Empty input vectors, no force will be computed" << "\n";
        return false;
    }
    
    if (x.size() != x0.size()) {
        msg_warning("BaseBeamHookeLawForceField") << "Position vector size (" << x.size() 
                                                << ") doesn't match rest position size (" << x0.size() << ")" << "\n";
        return false;
    }
    
    unsigned int sz = d_length.getValue().size();
    if (x.size() != sz) {
        msg_warning("BaseBeamHookeLawForceField") << "Length vector size (" << sz 
                                                << ") should have the same size as position vector (" << x.size() << ")" << "\n";
        return false;
    }
    
    if (d_variantSections.getValue()) {
        if (d_youngModulusList.getValue().size() != x.size()) {
            msg_warning("BaseBeamHookeLawForceField") << "Using variant sections but youngModulusList size (" 
                                                   << d_youngModulusList.getValue().size() 
                                                   << ") doesn't match position vector size (" << x.size() << ")" << "\n";
            return false;
        }
        
        if (d_poissonRatioList.getValue().size() != x.size()) {
            msg_warning("BaseBeamHookeLawForceField") << "Using variant sections but poissonRatioList size (" 
                                                   << d_poissonRatioList.getValue().size() 
                                                   << ") doesn't match position vector size (" << x.size() << ")" << "\n";
            return false;
        }
        
        if (m_K_sectionList.size() != x.size()) {
            msg_warning("BaseBeamHookeLawForceField") << "Using variant sections but section list size (" 
                                                   << m_K_sectionList.size() 
                                                   << ") doesn't match position vector size (" << x.size() << ")" << "\n";
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Validates material parameters
 * 
 * Ensures that Young's modulus and Poisson ratio have physically valid values.
 * 
 * @return true if parameters are valid, false otherwise
 */
template<typename DataTypes>
bool BaseBeamHookeLawForceField<DataTypes>::validateMaterialParameters() const
{
    // Check Young's modulus (must be positive)
    if (d_youngModulus.getValue() <= 0.0) {
        msg_error("BaseBeamHookeLawForceField") << "Young's modulus must be positive: " << d_youngModulus.getValue();
        return false;
    }
    
    // Check Poisson ratio (theoretical limits: -1.0 < nu < 0.5)
    if (d_poissonRatio.getValue() <= -1.0 || d_poissonRatio.getValue() >= 0.5) {
        msg_error("BaseBeamHookeLawForceField") << "Poisson ratio must be in range (-1.0, 0.5): " << d_poissonRatio.getValue();
        return false;
    }
    
    // For variant sections, check all values
    if (d_variantSections.getValue()) {
        const Vector& youngModulusList = d_youngModulusList.getValue();
        const Vector& poissonRatioList = d_poissonRatioList.getValue();
        
        for (size_t i = 0; i < youngModulusList.size(); ++i) {
            if (youngModulusList[i] <= 0.0) {
                msg_error("BaseBeamHookeLawForceField") << "Young's modulus in list at index " << i << " must be positive: " << youngModulusList[i];
                return false;
            }
        }
        
        for (size_t i = 0; i < poissonRatioList.size(); ++i) {
            if (poissonRatioList[i] <= -1.0 || poissonRatioList[i] >= 0.5) {
                msg_error("BaseBeamHookeLawForceField") << "Poisson ratio in list at index " << i << " must be in range (-1.0, 0.5): " << poissonRatioList[i];
                return false;
            }
        }
    }
    
    // If using inertia parameters directly, validate them
    if (d_useInertiaParams.getValue()) {
        // Only check these if they're actually provided (not default initialized)
        if (d_EIy.isSet() && d_EIy.getValue() <= 0.0) {
            msg_error("BaseBeamHookeLawForceField") << "EIy parameter must be positive: " << d_EIy.getValue();
            return false;
        }
        
        if (d_EIz.isSet() && d_EIz.getValue() <= 0.0) {
            msg_error("BaseBeamHookeLawForceField") << "EIz parameter must be positive: " << d_EIz.getValue();
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Validates cross-section parameters
 * 
 * Ensures that cross-section dimensions are physically valid.
 * 
 * @return true if parameters are valid, false otherwise
 */
template<typename DataTypes>
bool BaseBeamHookeLawForceField<DataTypes>::validateSectionParameters() const
{
    // Get cross-section shape
    const std::string& shape = d_crossSectionShape.getValue().getSelectedItem();
    
    if (shape == "circular") {
        // Check radius (must be positive)
        if (d_radius.getValue() <= 0.0) {
            msg_error("BaseBeamHookeLawForceField") << "External radius must be positive: " << d_radius.getValue();
            return false;
        }
        
        // Check inner radius (must be non-negative and less than external radius)
        if (d_innerRadius.getValue() < 0.0) {
            msg_error("BaseBeamHookeLawForceField") << "Inner radius cannot be negative: " << d_innerRadius.getValue();
            return false;
        }
        
        if (d_innerRadius.getValue() >= d_radius.getValue()) {
            msg_error("BaseBeamHookeLawForceField") << "Inner radius (" << d_innerRadius.getValue() 
                << ") must be less than external radius (" << d_radius.getValue() << ")";
            return false;
        }
    } 
    else if (shape == "rectangular") {
        // Check rectangular dimensions (must be positive)
        if (d_lengthY.getValue() <= 0.0) {
            msg_error("BaseBeamHookeLawForceField") << "Side length Y must be positive: " << d_lengthY.getValue();
            return false;
        }
        
        if (d_lengthZ.getValue() <= 0.0) {
            msg_error("BaseBeamHookeLawForceField") << "Side length Z must be positive: " << d_lengthZ.getValue();
            return false;
        }
    }
    else {
        msg_error("BaseBeamHookeLawForceField") << "Unknown cross-section shape: " << shape;
        return false;
    }
    
    // Check section lengths
    const Vector& lengths = d_length.getValue();
    for (size_t i = 0; i < lengths.size(); ++i) {
        if (lengths[i] <= 0.0) {
            msg_error("BaseBeamHookeLawForceField") << "Section length at index " << i << " must be positive: " << lengths[i];
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Computes cross-section properties
 * 
 * Calculates moment of inertia, cross-section area, and stiffness matrices
 * based on the cross-section shape and material properties.
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::computeSectionProperties()
{
    // Initialize matrices to zero
    m_K_section66.fill(0);
    m_K_sectionList.clear();
    m_K_section66List.clear();
    
    // Calculate moments of inertia based on cross-section shape
    if (d_crossSectionShape.getValue().getSelectedItem() == "rectangular") {
        // Rectangular cross-section
        const Real Ly = d_lengthY.getValue();
        const Real Lz = d_lengthZ.getValue();
        
        const Real LyLzLzLz = Ly * Lz * Lz * Lz;
        const Real LzLyLyLy = Lz * Ly * Ly * Ly;
        
        Iy = LyLzLzLz / 12.0;
        Iz = LzLyLyLy / 12.0;
        J = Iy + Iz;
        m_crossSectionArea = Ly * Lz;
    }
    else {
        // Circular cross-section
        const Real r = d_radius.getValue();
        const Real rInner = d_innerRadius.getValue();
        const Real r4 = r * r * r * r;
        const Real rInner4 = rInner * rInner * rInner * rInner;
        
        Iy = M_PI * (r4 - rInner4) / 4.0;
        Iz = Iy;
        J = Iy + Iz;
        m_crossSectionArea = M_PI * (r * r - rInner * rInner);
    }
    
    // Calculate stiffness matrices
    if (!d_variantSections.getValue()) {
        // Uniform section
        if (!d_useInertiaParams.getValue()) {
            // Calculate from material properties
            const Real E = d_youngModulus.getValue();
            const Real G = E / (2.0 * (1.0 + d_poissonRatio.getValue()));
            
            // 3x3 stiffness matrix for rotational DOFs
            m_K_section[0][0] = G * J;
            m_K_section[1][1] = E * Iy;
            m_K_section[2][2] = E * Iz;
            
            // 6x6 stiffness matrix if needed for more complex beam formulations
            if (DataTypes::spatial_dimensions > 3) {
                m_K_section66[0][0] = E * m_crossSectionArea;  // Axial stiffness
                m_K_section66[1][1] = G * m_crossSectionArea;  // Shear stiffness y
                m_K_section66[2][2] = G * m_crossSectionArea;  // Shear stiffness z
                m_K_section66[3][3] = G * J;                   // Torsional stiffness
                m_K_section66[4][4] = E * Iy;                 // Bending stiffness y
                m_K_section66[5][5] = E * Iz;                 // Bending stiffness z
            }
        }
        else {
            // Use user-provided inertia parameters
            msg_info("BaseBeamHookeLawForceField") << "Using pre-calculated inertia parameters for stiffness matrix.";
            
            // 3x3 stiffness matrix
            m_K_section[0][0] = d_GI.getValue();
            
            // Use specific EIy/EIz if provided, otherwise use general EI
            if (d_EIy.isSet()) {
                m_K_section[1][1] = d_EIy.getValue();
            } else {
                m_K_section[1][1] = d_EI.getValue();
            }
            
            if (d_EIz.isSet()) {
                m_K_section[2][2] = d_EIz.getValue();
            } else {
                m_K_section[2][2] = d_EI.getValue();
            }
            
            // 6x6 stiffness matrix if needed
            if (DataTypes::spatial_dimensions > 3) {
                m_K_section66[0][0] = d_EA.getValue();        // Axial stiffness
                m_K_section66[1][1] = d_GA.getValue();        // Shear stiffness y
                m_K_section66[2][2] = d_GA.getValue();        // Shear stiffness z
                m_K_section66[3][3] = d_GI.getValue();        // Torsional stiffness
                
                // Use specific EIy/EIz if provided, otherwise use general EI
                if (d_EIy.isSet()) {
                    m_K_section66[4][4] = d_EIy.getValue();
                } else {
                    m_K_section66[4][4] = d_EI.getValue();
                }
                
                if (d_EIz.isSet()) {
                    m_K_section66[5][5] = d_EIz.getValue();
                } else {
                    m_K_section66[5][5] = d_EI.getValue();
                }
            }
        }
    }
    else {
        // Variant sections
        msg_info("BaseBeamHookeLawForceField") << "Computing properties for variant beam sections.";
        
        const auto szL = d_length.getValue().size();
        if ((szL != d_poissonRatioList.getValue().size()) || (szL != d_youngModulusList.getValue().size())) {
            msg_error("BaseBeamHookeLawForceField") << "The size of data 'length', 'youngModulusList', and 'poissonRatioList' should be the same!";
            return;
        }
        
        for (auto k = 0; k < szL; k++) {
            // 3x3 stiffness matrix for each section
            Mat33 _m_K_section;
            const Real E = d_youngModulusList.getValue()[k];
            const Real G = E / (2.0 * (1.0 + d_poissonRatioList.getValue()[k]));
            
            _m_K_section[0][0] = G * J;
            _m_K_section[1][1] = E * Iy;
            _m_K_section[2][2] = E * Iz;
            m_K_sectionList.push_back(_m_K_section);
            
            // 6x6 stiffness matrix if needed
            if (DataTypes::spatial_dimensions > 3) {
                Mat66 _m_K_section66;
                _m_K_section66[0][0] = E * m_crossSectionArea;  // Axial stiffness
                _m_K_section66[1][1] = G * m_crossSectionArea;  // Shear stiffness y
                _m_K_section66[2][2] = G * m_crossSectionArea;  // Shear stiffness z
                _m_K_section66[3][3] = G * J;                   // Torsional stiffness
                _m_K_section66[4][4] = E * Iy;                 // Bending stiffness y
                _m_K_section66[5][5] = E * Iz;                 // Bending stiffness z
                m_K_section66List.push_back(_m_K_section66);
            }
        }
        
        if (d_useInertiaParams.getValue()) {
            msg_warning("BaseBeamHookeLawForceField") << "Using variant sections with pre-calculated inertia parameters is not yet supported.";
        }
    }
}

/******************************************************************************
* FORCE COMPUTATION METHODS                                                  *
******************************************************************************/

/**
 * @brief Adds forces for uniform section beams
 * 
 * Computes and adds forces to the force vector for beams with uniform cross-section
 * properties. Supports both single-threaded and multi-threaded computation.
 * 
 * @param d_f Force vector
 * @param d_x Position vector
 * @param d_x0 Rest position vector
 * @param lengths Vector of beam section lengths
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::addForceUniformSection(DataVecDeriv& d_f, 
                                                                  const DataVecCoord& d_x, 
                                                                  const DataVecCoord& d_x0,
                                                                  const type::vector<Real>& lengths)
{
    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    const VecCoord& x0 = d_x0.getValue();
    const size_t size = x.size();
    sofa::helper::ScopedAdvancedTimer timer("UniformSection");
    
    if (d_useMultiThreading.getValue() && size > 1) {
        sofa::simulation::TaskScheduler::Task::Status status;
        sofa::simulation::TaskScheduler& scheduler = *(sofa::simulation::TaskScheduler::getInstance());
        // Define a lambda for parallel execution
        auto calcForce = [this, &f, &x, &x0, &lengths](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                f[i] -= (m_K_section * (x[i] - x0[i])) * lengths[i];
            }
        };
        
        // Determine chunk size for parallel processing
        const size_t chunk_size = std::max<size_t>(1, size / (scheduler.getThreadCount() * 2));
        size_t start = 0;
        
        // Create and queue tasks for parallel execution
        while (start < size) {
            size_t end = std::min(start + chunk_size, size);
            Task* task = new Task(std::bind(calcForce, start, end), &status);
            scheduler.addTask(task);
            start = end;
        }
        
        // Wait for all tasks to complete
        scheduler.workUntilDone(&status);
    } 
    else {
        // Single-threaded fallback
        for (size_t i = 0; i < size; ++i) {
            f[i] -= (m_K_section * (x[i] - x0[i])) * lengths[i];
        }
    }
    
    d_f.endEdit();
}
/**
 * @brief Adds forces for variant section beams
 * 
 * Computes and adds forces to the force vector for beams with varying cross-section
 * properties. Supports both single-threaded and multi-threaded computation.
 * 
 * @param d_f Force vector
 * @param d_x Position vector
 * @param d_x0 Rest position vector
 * @param lengths Vector of beam section lengths
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::addForceVariantSection(DataVecDeriv& d_f, 
                                                                  const DataVecCoord& d_x, 
                                                                  const DataVecCoord& d_x0,
                                                                  const type::vector<Real>& lengths)
{
    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();
    const VecCoord& x0 = d_x0.getValue();
    const size_t size = x.size();
    sofa::helper::ScopedAdvancedTimer timer("VariantSection");
    
    if (d_useMultiThreading.getValue() && size > 1) {
        sofa::simulation::TaskScheduler::Task::Status status;
        sofa::simulation::TaskScheduler& scheduler = *(sofa::simulation::TaskScheduler::getInstance());
        
        // Define a lambda for parallel execution
        auto calcForce = [this, &f, &x, &x0, &lengths](size_t start, size_t end) {
            for (size_t i = start; i < end; ++i) {
                f[i] -= (m_K_sectionList[i] * (x[i] - x0[i])) * lengths[i];
            }
        };
        
        // Determine chunk size for parallel processing
        const size_t chunk_size = std::max<size_t>(1, size / (scheduler.getThreadCount() * 2));
        size_t start = 0;
        
        // Create and queue tasks for parallel execution
        while (start < size) {
            size_t end = std::min(start + chunk_size, size);
            Task* task = new Task(std::bind(calcForce, start, end), &status);
            scheduler.addTask(task);
            start = end;
        }
        
        // Wait for all tasks to complete
        scheduler.workUntilDone(&status);
    } 
    else {
        // Single-threaded fallback
        for (size_t i = 0; i < size; ++i) {
            f[i] -= (m_K_sectionList[i] * (x[i] - x0[i])) * lengths[i];
        }
    }
    
    d_f.endEdit();
}
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* mparams,
                                                   DataVecDeriv& d_f,
                                                   const DataVecCoord& d_x,
                                                   const DataVecDeriv& d_v)
{
    SOFA_UNUSED(d_v);
    SOFA_UNUSED(mparams);
    
    sofa::helper::ScopedAdvancedTimer timer("BaseBeamHookeLawForceField::addForce");
    
    // Get rest position
    const DataVecCoord& d_x0 = *this->mstate->read(sofa::core::vec_id::read_access::restPosition);
    
    // Validate input data
    if (!validateInputData(d_f, d_x, d_x0)) {
        compute_df = false;
        return;
    }
    
    // Prepare force vector
    VecDeriv& f = *d_f.beginEdit();
    f.resize(d_x.getValue().size());
    d_f.endEdit();
    
    const type::vector<Real>& lengths = d_length.getValue();
    
    // Call the appropriate specialized force calculation method
    if (!d_variantSections.getValue()) {
        addForceUniformSection(d_f, d_x, d_x0, lengths);
    } else {
        addForceVariantSection(d_f, d_x, d_x0, lengths);
    }
}

/**
 * @brief Computes the product of the stiffness matrix and a displacement vector
 * 
 * This method is used for implicit integration schemes to compute the change in force
 * due to a displacement.
 * 
 * @param mparams Mechanical parameters
 * @param d_df Differential force vector
 * @param d_dx Differential displacement vector
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams,
                                                    DataVecDeriv& d_df,
                                                    const DataVecDeriv& d_dx)
{
    if (!compute_df)
        return;

    WriteAccessor<DataVecDeriv> df = d_df;
    ReadAccessor<DataVecDeriv> dx = d_dx;
    const Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    df.resize(dx.size());
    if (!d_variantSections.getValue()) {
        // For uniform section beams
        for (unsigned int i = 0; i < dx.size(); i++) {
            df[i] -= (m_K_section * dx[i]) * kFactor * d_length.getValue()[i];
        }
    }
    else {
        // For variant section beams
        for (unsigned int i = 0; i < dx.size(); i++) {
            df[i] -= (m_K_sectionList[i] * dx[i]) * kFactor * d_length.getValue()[i];
        }
    }
}

/******************************************************************************
* MATRIX COMPUTATION METHODS                                                 *
******************************************************************************/

/**
 * @brief Adds the stiffness matrix contribution to the global matrix
 * 
 * This method is called during the assembly of the global stiffness matrix to add
 * the contribution of this force field.
 * 
 * @param mparams Mechanical parameters
 * @param matrix Multi-matrix accessor
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::addKToMatrix(const sofa::core::MechanicalParams* mparams,
                                                       const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    const VecCoord& pos = this->mstate->read(sofa::core::VecCoordId::position())->getValue();
    addKToMatrixImpl(mat, offset, pos, d_length.getValue(), kFact);
}

/**
 * @brief Template implementation to add stiffness to matrix
 * 
 * Implementation details for adding stiffness to the global matrix.
 * This is specialized in derived classes for different data types.
 * 
 * @param mat Base matrix
 * @param offset Matrix offset
 * @param pos Position vector
 * @param lengths Vector of beam section lengths
 * @param kFact Stiffness factor
 */
template<typename DataTypes>
template<class MatrixType>
void BaseBeamHookeLawForceField<DataTypes>::addKToMatrixImpl(MatrixType* mat, 
                                                           unsigned int offset, 
                                                           const VecCoord& pos, 
                                                           const Vector& lengths, 
                                                           Real kFact)
{
    for (unsigned int n = 0; n < pos.size(); n++)
    {
        if (!d_variantSections.getValue()) {
            // For uniform section beams
            for (unsigned int i = 0; i < 3; i++) {
                for (unsigned int j = 0; j < 3; j++) {
                    mat->add(offset + i + 3*n, offset + j + 3*n, -kFact * m_K_section[i][j] * lengths[n]);
                }
            }
        }
        else {
            // For variant section beams
            for (unsigned int i = 0; i < 3; i++) {
                for (unsigned int j = 0; j < 3; j++) {
                    mat->add(offset + i + 3*n, offset + j + 3*n, -kFact * m_K_sectionList[n][i][j] * lengths[n]);
                }
            }
        }
    }
}

/******************************************************************************
* UTILITY METHODS                                                            *
******************************************************************************/

/**
 * @brief Computes the potential energy of the beam
 * 
 * This method returns the strain energy stored in the beam due to deformation.
 * The strain energy is calculated as 0.5 * displacement^T * K * displacement.
 * 
 * @param mparams Mechanical parameters
 * @param d_x Position vector
 * @return Potential energy
 */
template<typename DataTypes>
double BaseBeamHookeLawForceField<DataTypes>::getPotentialEnergy(const sofa::core::MechanicalParams* mparams,
                                                               const DataVecCoord& d_x) const
{
    SOFA_UNUSED(mparams);
    
    double energy = 0.0;
    const VecCoord& x = d_x.getValue();
    const VecCoord& x0 = this->mstate->read(sofa::core::VecCoordId::restPosition())->getValue();
    const type::vector<Real>& lengths = d_length.getValue();
    
    // Make sure we have valid data
    if (x.size() != x0.size() || x.size() != lengths.size()) {
        return 0.0;
    }
    
    // Calculate the potential energy
    if (!d_variantSections.getValue()) {
        // For uniform section beams
        for (unsigned int i = 0; i < x.size(); i++) {
            Deriv delta = x[i] - x0[i];
            energy += 0.5 * (delta * (m_K_section * delta)) * lengths[i];
        }
    }
    else {
        // For variant section beams
        for (unsigned int i = 0; i < x.size(); i++) {
            Deriv delta = x[i] - x0[i];
            energy += 0.5 * (delta * (m_K_sectionList[i] * delta)) * lengths[i];
        }
    }
    
    return energy;
}

/**
 * @brief Returns the external radius of the beam
 * 
 * Utility method to access the beam's external radius.
 * 
 * @return The external radius of the beam
 */
template<typename DataTypes>
typename BaseBeamHookeLawForceField<DataTypes>::Real BaseBeamHookeLawForceField<DataTypes>::getRadius()
{
    return d_radius.getValue();
}

/**
 * @brief Compute relative rotations between beam cross-sections
 * 
 * Updates m_relativeRotations based on the current frames.
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::computeRelativeRotations()
{
    const size_t numSegments = m_currentFrames.size() - 1;
    
    for (size_t i = 0; i < numSegments; ++i) {
        // Compute the relative rotation from frame i to frame i+1
        m_relativeRotations[i] = m_currentFrames[i].rotation().inverse() * 
                                m_currentFrames[i+1].rotation();
    }
}

/**
 * @brief Update current frames based on positions and rotations
 * 
 * @param positions Current positions of the beam nodes
 * @param rotations Current rotations of the beam cross-sections
 */
template<typename DataTypes>
void BaseBeamHookeLawForceField<DataTypes>::updateFrames(
    const VecCoord& positions, 
    const std::vector<SO3Type>& rotations)
{
    const size_t numNodes = positions.size();
    assert(rotations.size() == numNodes);
    
    for (size_t i = 0; i < numNodes; ++i) {
        // Extract position
        Vector3 position = getPosition(positions[i]);
        
        // Update current frame
        m_currentFrames[i] = SE3Type(rotations[i], position);
    }
    
    // Update relative rotations
    computeRelativeRotations();
}

/**
 * @brief Get local frame axis vectors
 * @param frameIndex Index of the frame
 * @return Tuple of local x, y, z axes as Vector3
 */
template<typename DataTypes>
std::tuple<typename BaseBeamHookeLawForceField<DataTypes>::Vector3,
           typename BaseBeamHookeLawForceField<DataTypes>::Vector3,
           typename BaseBeamHookeLawForceField<DataTypes>::Vector3>
BaseBeamHookeLawForceField<DataTypes>::getLocalAxes(size_t frameIndex) const
{
    using Matrix3 = typename SO3Type::Matrix;
    
    // Get rotation matrix
    const Matrix3 R = m_currentFrames[frameIndex].rotation().matrix();
    
    // Extract column vectors
    Vector3 xAxis(R(0,0), R(1,0), R(2,0));
    Vector3 yAxis(R(0,1), R(1,1), R(2,1));
    Vector3 zAxis(R(0,2), R(1,2), R(2,2));
    
    return std::make_tuple(xAxis, yAxis, zAxis);
}

/**
 * @brief Compute the strain tensor between two frames
 * @param frame1 First frame
 * @param frame2 Second frame
 * @return Strain tensor in se(3)
 */
template<typename DataTypes>
typename BaseBeamHookeLawForceField<DataTypes>::SE3TangentType
BaseBeamHookeLawForceField<DataTypes>::computeStrainTensor(
    const SE3Type& frame1,
    const SE3Type& frame2) const
{
    // Compute relative transformation
    SE3Type relativeTransform = frame1.inverse() * frame2;
    
    // Extract strain using logarithm mapping to the Lie algebra
    return relativeTransform.log();
}

/**
 * @brief Transform a local force/moment to the global frame
 * 
 * @param localWrench Local force/moment in the cross-section frame
 * @param frameIndex Index of the cross-section frame
 * @return Global force/moment
 */
template<typename DataTypes>
typename BaseBeamHookeLawForceField<DataTypes>::Deriv 
BaseBeamHookeLawForceField<DataTypes>::transformWrenchToGlobal(
    const Deriv& localWrench, 
    size_t frameIndex) const
{
    // Extract force and moment components
    Vector3 localForce = getForce(localWrench);
    Vector3 localMoment = getMoment(localWrench);
    
    // Transform to global frame
    Vector3 globalForce = m_currentFrames[frameIndex].rotation().act(localForce);
    Vector3 globalMoment = m_currentFrames[frameIndex].rotation().act(localMoment);
    
    // Reconstruct Deriv
    return createDeriv(globalForce, globalMoment);
}

/**
 * @brief Transform a global force/moment to a local frame
 * 
 * @param globalWrench Global force/moment
 * @param frameIndex Index of the cross-section frame
 * @return Local force/moment
 */
template<typename DataTypes>
typename BaseBeamHookeLawForceField<DataTypes>::Deriv 
BaseBeamHookeLawForceField<DataTypes>::transformWrenchToLocal(
    const Deriv& globalWrench, 
    size_t frameIndex) const
{
    // Extract force and moment components
    Vector3 globalForce = getForce(globalWrench);
    Vector3 globalMoment = getMoment(globalWrench);
    
    // Transform to local frame
    Vector3 localForce = m_currentFrames[frameIndex].rotation().inverse().act(globalForce);
    Vector3 localMoment = m_currentFrames[frameIndex].rotation().inverse().act(globalMoment);
    
    // Reconstruct Deriv
    return createDeriv(localForce, localMoment);
}

/**
 * @brief Extract strains using Lie Group formalism
 * 
 * Computes the strains (curvature/twist and stretch) between consecutive frames
 * using the Lie algebra logarithm.
 * 
 * @return Vector of strains for each beam segment
 */
template<typename DataTypes>
std::vector<typename BaseBeamHookeLawForceField<DataTypes>::Deriv> 
BaseBeamHookeLawForceField<DataTypes>::computeStrains() const
{
    const size_t numSegments = m_currentFrames.size() - 1;
    std::vector<Deriv> strains(numSegments);
    
    for (size_t i = 0; i < numSegments; ++i) {
        // Compute relative transformation from reference to current configuration
        SE3Type refRelTransform = m_referenceFrames[i].inverse() * m_referenceFrames[i+1];
        SE3Type curRelTransform = m_currentFrames[i].inverse() * m_currentFrames[i+1];
        
        // Compute deformation (reference to current)
        SE3Type deformation = refRelTransform.inverse() * curRelTransform;
        
        // Extract strain from the Lie algebra using logarithm
        SE3TangentType se3Strain = deformation.log();
        
        // Convert to beam strain representation
        // First 3 components are translation strain, last 3 are rotation strain
        Vector3 transStrain(se3Strain[0], se3Strain[1], se3Strain[2]);
        Vector3 rotStrain(se3Strain[3], se3Strain[4], se3Strain[5]);
        
        // Scale by segment length
        const Real segmentLength = d_length.getValue()[i];
        transStrain /= segmentLength;
        rotStrain /= segmentLength;
        
        // Create strain Deriv
        strains[i] = createDeriv(transStrain, rotStrain);
    }
    
    return strains;
}

/**
 * @brief Get the current frame transformation for a specific cross-section
 * 
 * @param index Cross-section index
 * @
