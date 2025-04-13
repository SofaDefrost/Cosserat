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
#include <Cosserat/forcefield/BeamHookeLawForceField.h>

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
 * @brief Implementation of all methods for the BeamHookeLawForceField class
 * 
 * This force field simulates the mechanical behavior of beam elements using Hooke's law.
 * It supports both uniform sections (same properties throughout the beam) and variant sections
 * (different properties for different segments).
 */

/******************************************************************************
* CONSTRUCTORS / DESTRUCTORS                                                 *
******************************************************************************/

template<typename DataTypes>
BeamHookeLawForceField<DataTypes>::BeamHookeLawForceField()
    : Inherit1(),
    d_crossSectionShape( initData(&d_crossSectionShape, {"circular","rectangular"},
                               "crossSectionShape",
                               "shape of the cross-section. Can be: circular (tube with external radius being radius and internal radius being innerRadius ) or rectangular (lengthY and lengthZ) . Default is circular" )),
    d_youngModulus( initData( &d_youngModulus, 1.0e9, "youngModulus", "Young Modulus describes the stiffness of the material")),
    d_poissonRatio( initData( &d_poissonRatio, 0.45, "poissonRatio", "poisson Ratio describes the compressibility of the material")),
    d_length( initData( &d_length, "length", "The list of lengths of the different beam's sections.")),
    d_radius( initData( &d_radius, 1.0, "radius", "external radius of the cross section (if circular)")),
    d_innerRadius( initData( &d_innerRadius, 0.0, "innerRadius", "internal radius of the cross section (if circular)")),
    d_lengthY( initData( &d_lengthY, 1.0, "lengthY", "side length of the cross section along local y axis (if rectangular)")),
    d_lengthZ( initData( &d_lengthZ, 1.0, "lengthZ", "side length of the cross section along local z axis (if rectangular)")),
    d_variantSections(initData(&d_variantSections, false, "variantSections", "In case we have variant beam sections this has to be set to true")),
    d_youngModulusList(initData(&d_youngModulusList, "youngModulusList", "The list of Young modulus in case we have sections with variable physical properties")),
    d_poissonRatioList(initData(&d_poissonRatioList, "poissonRatioList", "The list of poisson's ratio in case we have sections with variable physical properties")),
    d_useInertiaParams(initData(&d_useInertiaParams, false, "useInertiaParams", "If the inertia parameters are given by the user, there is no longer any need to use @d_youngModulus and @d_poissonRatio.")),
    d_GI(initData(&d_GI, "GI", "The inertia parameter, GI")),
    d_GA(initData(&d_GA, "GA", "The inertia parameter, GA")),
    d_EA(initData(&d_EA, "EA", "The inertia parameter, EA")),
    d_EI(initData(&d_EI, "EI", "The inertia parameter, EI")),
    d_useMultiThreading(initData(&d_useMultiThreading, false, "useMultiThreading", "Use multithreading for force computation"))
{
    compute_df = true;
}

template<typename DataTypes>
BeamHookeLawForceField<DataTypes>::~BeamHookeLawForceField() = default;

/******************************************************************************
* INITIALIZATION METHODS                                                     *
******************************************************************************/

/**
 * @brief Initialize the force field, setting up parameters and cross-section properties
 * 
 * This method is called during initialization to set up the force field.
 * It delegates to reinit() to calculate cross-section properties.
 */
template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::init()
{
    Inherit1::init();
    reinit();
}

/**
 * @brief Reinitialize the force field, recalculating cross-section properties
 * 
 * This method calculates the moments of inertia and stiffness matrices
 * based on the chosen cross-section shape and material properties.
 * It supports both uniform sections and variant sections.
 */
template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::reinit()
{
    // Precompute and store inertia values
    Real Iy, Iz, J, A;
    if (d_crossSectionShape.getValue().getSelectedItem() == "rectangular")  // rectangular cross-section
    {
        Real Ly = d_lengthY.getValue();
        Real Lz = d_lengthZ.getValue();

        const Real LyLzLzLz = Ly * Lz * Lz * Lz;
        const Real LzLyLyLy = Lz * Ly * Ly * Ly;

        Iy = LyLzLzLz / 12.0;
        Iz = LzLyLyLy / 12.0;
        J = Iy + Iz;
        A = Ly * Lz;

    }
    else // circular cross-section
    {
        msg_info() << "Cross section shape." << d_crossSectionShape.getValue().getSelectedItem();

        Real r = d_radius.getValue();
        Real rInner = d_innerRadius.getValue();
        const Real r4 = r * r * r * r;
        const Real rInner4 = rInner * rInner * rInner * rInner;

        Iy = M_PI * (r4 - rInner4) / 4.0;
        Iz = Iy;
        J = Iy + Iz;
        A = M_PI * (r * r - rInner * rInner);

    }
    m_crossSectionArea = A;

    // if we are dealing with different physical properties: YM and PR
    if (!d_variantSections.getValue())
    {
        if (!d_useInertiaParams.getValue())
        {
            Real E = d_youngModulus.getValue();
            Real G = E / (2.0 * (1.0 + d_poissonRatio.getValue()));
            // Inertial matrix
            m_K_section[0][0] = G * J;
            m_K_section[1][1] = E * Iy;
            m_K_section[2][2] = E * Iz;
        }
        else
        {
            msg_info("BeamHookeLawForceField") << "Pre-calculated inertia parameters are used for the computation "
                                              "of the stiffness matrix.";
            m_K_section[0][0] = d_GI.getValue();
            m_K_section[1][1] = d_EI.getValue();
            m_K_section[2][2] = d_EI.getValue();
        }

    }
    else {
        msg_info("BeamHookeLawForceField") << "Multi section beam are used for the simulation!";
        m_K_sectionList.clear();

        const auto szL = d_length.getValue().size();
        if ((szL != d_poissonRatioList.getValue().size()) || (szL != d_youngModulusList.getValue().size())) {
            msg_error("BeamHookeLawForceField") << "Please the size of the data length, youngModulusList and "
                                               "poissonRatioList should be the same!";
            return;
        }

        for (auto k = 0; k < szL; k++)
        {
            Mat33 _m_K_section;
            Real E = d_youngModulusList.getValue()[k];
            Real G = E / (2.0 * (1.0 + d_poissonRatioList.getValue()[k]));

            _m_K_section[0][0] = G * J;
            _m_K_section[1][1] = E * Iy;
            _m_K_section[2][2] = E * Iz;
            m_K_sectionList.push_back(_m_K_section);
        }
        msg_info("BeamHookeLawForceField") << "If you plan to use a multi section beam with (different "
                                          "mechanical properties) and pre-calculated inertia parameters "
                                          "(GI, GA, etc.), this is not yet supported.";
    }
}

/**
 * @brief Validates input data for force computation
 * 
 * Checks that the input data is valid for force computation, including:
 * - The mechanical state exists
 * - Position vectors are not empty and have matching sizes
 * - The length vector has the correct size
 * - For variant sections, the section list has the correct size
 * 
 * @param d_f Force vector
 * @param d_x Position vector
 * @param d_x0 Rest position vector
 * @return true if data is valid, false otherwise
 */
template<typename DataTypes>
bool BeamHookeLawForceField<DataTypes>::validateInputData(const DataVecDeriv& d_f, 
                                                          const DataVecCoord& d_x, 
                                                          const DataVecCoord& d_x0) const
{
    if (!this->getMState()) {
        msg_info("BeamHookeLawForceField") << "No Mechanical State found, no force will be computed..." << "\n";
        return false;
    }
    
    const VecCoord& x = d_x.getValue();
    const VecCoord& x0 = d_x0.getValue();
    
    if (x.empty() || x0.empty()) {
        msg_warning("BeamHookeLawForceField") << "Empty input vectors, no force will be computed" << "\n";
        return false;
    }
    
    if (x.size() != x0.size()) {
        msg_warning("BeamHookeLawForceField") << "Position vector size (" << x.size() 
                                             << ") doesn't match rest position size (" << x0.size() << ")" << "\n";
        return false;
    }
    
    unsigned int sz = d_length.getValue().size();
    if (x.size() != sz) {
        msg_warning("BeamHookeLawForceField") << "Length vector size (" << sz 
                                             << ") should have the same size as position vector (" << x.size() << ")" << "\n";
        return false;
    }
    
    if (d_variantSections.getValue() && m_K_sectionList.size() != x.size()) {
        msg_warning("BeamHookeLawForceField") << "Using variant sections but section list size (" << m_K_sectionList.size() 
                                             << ") doesn't match position vector size (" << x.size() << ")" << "\n";
        return false;
    }
    
    return true;
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
void BeamHookeLawForceField<DataTypes>::addForceUniformSection(DataVecDeriv& d_f, 
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
void BeamHookeLawForceField<DataTypes>::addForceVariantSection(DataVecDeriv& d_f, 
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
template<typename DataTypes>
void BeamHookeLawForceField<DataTypes>::addForce(const MechanicalParams* mparams,
                                                 DataVecDeriv& d_f,
                                                 const DataVecCoord& d_x,
                                                 const DataVecDeriv& d_v)
{
    SOFA_UNUSED(d_v);
    SOFA_UNUSED(mparams);
    
    sofa::helper::ScopedAdvancedTimer timer("BeamHookeLawForceField::addForce");
    
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
void BeamHookeLawForceField<DataTypes>::addDForce(const MechanicalParams* mparams,
                                                  DataVecDeriv& d_df,
                                                  const DataVecDeriv& d_dx)
{
    if (!compute_df)
        return;

    WriteAccessor< DataVecDeriv > df = d_df;
    ReadAccessor< DataVecDeriv > dx = d_dx;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    df.resize(dx.size());
    if (!d_variantSections.getValue()) {
        for (unsigned int i = 0; i < dx.size(); i++) {
            df[i] -= (m_K_section * dx[i]) * kFactor * d_length.getValue()[i];
        }
    }
    else {
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
void BeamHookeLawForceField<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                     const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const VecCoord& pos = this->mstate->read(core::vec_id::read_access::position)->getValue();
    for (unsigned int n = 0; n < pos.size(); n++)
    {
        if (!d_variantSections.getValue()) {
            for (unsigned int i = 0; i < 3; i++) {
                for (unsigned int j = 0; j < 3; j++) {
                    mat->add(offset + i + 3*n, offset + j + 3*n, -kFact * m_K_section[i][j] * d_length.getValue()[n]);
                }
            }
        }
        else {
            for (unsigned int i = 0; i < 3; i++) {
                for (unsigned int j = 0; j < 3; j++) {
                    mat->add(offset + i + 3*n, offset + j + 3*n, -kFact * m_K_sectionList[n][i][j] * d_length.getValue()[n]);
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
 * Currently returns 0 as it is not implemented.
 * 
 * @param mparams Mechanical parameters
 * @param d_x Position vector
 * @return Potential energy
 */
template<typename DataTypes>
double BeamHookeLawForceField<DataTypes>::getPotentialEnergy(const MechanicalParams* mparams,
                                                             const DataVecCoord& d_x) const
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(d_x);

    return 0.0;
}

/**
 * @brief Returns the external radius of the beam
 * 
 * Utility method to access the beam's external radius.
 * 
 * @return The external radius of the beam
 */
template<typename DataTypes>
typename BeamHookeLawForceField<DataTypes>::Real BeamHookeLawForceField<DataTypes>::getRadius()
{
    return d_radius.getValue();
}

} // namespace sofa::component::forcefield
