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

#include <sofa/core/behavior/ForceField.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>
#include <sofa/simulation/TaskScheduler.h>

namespace sofa::component::forcefield
{

/**
 * @brief Base class for beam force fields using Hooke's law
 * 
 * This force field simulates the mechanical behavior of beam elements using Hooke's law.
 * It provides core functionality for both uniform sections (same properties throughout the beam)
 * and variant sections (different properties for different segments).
 * 
 * The class handles both circular and rectangular cross-sections, and supports
 * computation of beam properties either from material properties (Young's modulus,
 * Poisson ratio) or from directly specified inertia parameters.
 */
template<class DataTypes>
class BaseBeamHookeLawForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BaseBeamHookeLawForceField, DataTypes), SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));
    
    typedef sofa::core::behavior::ForceField<DataTypes> Inherit1;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef Data<VecCoord> DataVecCoord;
    typedef Data<VecDeriv> DataVecDeriv;
    typedef Data<MatrixDeriv> DataMatrixDeriv;
    typedef sofa::core::behavior::MultiMatrixAccessor MultiMatrixAccessor;
    typedef sofa::linearalgebra::BaseMatrix BaseMatrix;
    typedef sofa::type::vector<Real> Vector;
    
    // For task scheduling in multi-threaded computation
    typedef sofa::simulation::TaskScheduler::Task Task;
    
    // Matrix types for stiffness representation
    typedef sofa::type::Mat<3, 3, Real> Mat33;
    typedef sofa::type::Mat<6, 6, Real> Mat66;
    typedef sofa::type::vector<Mat33> VecMat33;
    typedef sofa::type::vector<Mat66> VecMat66;
    
protected:
    // Default constructor
    BaseBeamHookeLawForceField();
    
    // Destructor
    virtual ~BaseBeamHookeLawForceField();
    
    // Cross-section inertia values
    Real Iy, Iz, J;
    
    // Stiffness matrices
    Mat33 m_K_section;        // 3x3 stiffness matrix for rotational DOFs
    Mat66 m_K_section66;      // 6x6 stiffness matrix for rotational and translational DOFs
    VecMat33 m_K_sectionList; // List of stiffness matrices for variant sections
    VecMat66 m_K_section66List; // List of 6x6 stiffness matrices for variant sections
    
    // Flag to track whether to compute derivative forces
    bool compute_df;
    
    // Flag to track initialization state
    bool m_initialized;
    
    // Cross-section area
    Real m_crossSectionArea;
    
public:
    //////////////////////////////////////////////////////////////////////////
    // DATA MEMBERS
    //////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Shape of the cross-section
     * Can be "circular" (tube with external radius being radius and internal radius 
     * being innerRadius) or "rectangular" (lengthY and lengthZ)
     */
    Data<sofa::helper::OptionsGroup> d_crossSectionShape;
    
    /// Young's modulus - describes the stiffness of the material
    Data<Real> d_youngModulus;
    
    /// Poisson's ratio - describes the compressibility of the material
    Data<Real> d_poissonRatio;
    
    /// List of lengths of the different beam's sections
    Data<Vector> d_length;
    
    /// External radius of the cross section (if circular)
    Data<Real> d_radius;
    
    /// Internal radius of the cross section (if circular)
    Data<Real> d_innerRadius;
    
    /// Side length of the cross section along local y axis (if rectangular)
    Data<Real> d_lengthY;
    
    /// Side length of the cross section along local z axis (if rectangular)
    Data<Real> d_lengthZ;
    
    /// Set to true if we have variant beam sections
    Data<bool> d_variantSections;
    
    /// List of Young's modulus values for variant sections
    Data<Vector> d_youngModulusList;
    
    /// List of Poisson's ratio values for variant sections
    Data<Vector> d_poissonRatioList;
    
    /// Set to true to use provided inertia parameters instead of Young's modulus and Poisson ratio
    Data<bool> d_useInertiaParams;
    
    /// The inertia parameter GI (torsional rigidity)
    Data<Real> d_GI;
    
    /// The inertia parameter GA (shear rigidity)
    Data<Real> d_GA;
    
    /// The inertia parameter EA (axial rigidity)
    Data<Real> d_EA;
    
    /// The inertia parameter EI (bending rigidity)
    Data<Real> d_EI;
    
    /// The inertia parameter EIy (bending rigidity around y-axis)
    Data<Real> d_EIy;
    
    /// The inertia parameter EIz (bending rigidity around z-axis)
    Data<Real> d_EIz;
    
    /// Enable parallel processing for force computation
    Data<bool> d_useMultiThreading;
    
    //////////////////////////////////////////////////////////////////////////
    // METHODS
    //////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Initializes the force field
     * 
     * This method is called during initialization to set up the force field parameters.
     * It validates input data and computes cross-section properties.
     */
    virtual void init() override;
    
    /**
     * @brief Reinitializes the force field
     * 
     * Recalculates cross-section properties and stiffness matrices based on current parameters.
     */
    virtual void reinit() override;
    
    /**
     * @brief Validates input data for force computation
     * 
     * Checks that the input data is valid for force computation, including:
     * - The mechanical state exists
     * - Position vectors are not empty and have matching sizes
     * - The length vector has the correct size
     * - For variant sections, the section lists have the correct sizes
     * 
     * @param f Force vector
     * @param x Position vector
     * @param x0 Rest position vector
     * @return true if data is valid, false otherwise
     */
    virtual bool validateInputData(const DataVecDeriv& f, const DataVecCoord& x, const DataVecCoord& x0) const;
    
    /**
     * @brief Validates material parameters
     * 
     * Ensures that Young's modulus and Poisson ratio have physically valid values.
     * 
     * @return true if parameters are valid, false otherwise
     */
    virtual bool validateMaterialParameters() const;
    
    /**
     * @brief Validates cross-section parameters
     * 
     * Ensures that cross-section dimensions are physically valid.
     * 
     * @return true if parameters are valid, false otherwise
     */
    virtual bool validateSectionParameters() const;
    
    /**
     * @brief Computes cross-section properties
     * 
     * Calculates moment of inertia, cross-section area, and stiffness matrices
     * based on the cross-section shape and material properties.
     */
    virtual void computeSectionProperties();
    
    /**
     * @brief Adds forces for uniform section beams
     * 
     * Computes and adds forces to the force vector for beams with uniform cross-section
     * properties. Supports both single-threaded and multi-threaded computation.
     * 
     * @param f Force vector
     * @param x Position vector
     * @param x0 Rest position vector
     * @param lengths Vector of beam section lengths
     */
    virtual void addForceUniformSection(DataVecDeriv& f, const DataVecCoord& x, const DataVecCoord& x0, 
                                       const type::vector<Real>& lengths);
    
    /**
     * @brief Adds forces for variant section beams
     * 
     * Computes and adds forces to the force vector for beams with varying cross-section
     * properties. Supports both single-threaded and multi-threaded computation.
     * 
     * @param f Force vector
     * @param x Position vector
     * @param x0 Rest position vector
     * @param lengths Vector of beam section lengths
     */
    virtual void addForceVariantSection(DataVecDeriv& f, const DataVecCoord& x, const DataVecCoord& x0, 
                                       const type::vector<Real>& lengths);
    
    /**
     * @brief Main method for force computation
     * 
     * This method dispatches the force computation to either addForceUniformSection or
     * addForceVariantSection based on the beam configuration.
     * 
     * @param mparams Mechanical parameters
     * @param f Force vector
     * @param x Position vector
     * @param v Velocity vector
     */
    virtual void addForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& f, 
                        const DataVecCoord& x, const DataVecDeriv& v) override;
    
    /**
     * @brief Computes the product of the stiffness matrix and a displacement vector
     * 
     * This method is used for implicit integration schemes to compute the change in force
     * due to a displacement.
     * 
     * @param mparams Mechanical parameters
     * @param df Differential force vector
     * @param dx Differential displacement vector
     */
    virtual void addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& df, 
                         const DataVecDeriv& dx) override;
    
    /**
     * @brief Adds the stiffness matrix contribution to the global matrix
     * 
     * This method is called during the assembly of the global stiffness matrix to add
     * the contribution of this force field.
     * 
     * @param mparams Mechanical parameters
     * @param matrix Multi-matrix accessor
     */
    virtual void addKToMatrix(const sofa::core::MechanicalParams* mparams, 
                            const MultiMatrixAccessor* matrix) override;
    
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
    template<class MatrixType>
    void addKToMatrixImpl(MatrixType* mat, unsigned int offset, const VecCoord& pos, 
                         const Vector& lengths, Real kFact);
    
    /**
     * @brief Computes the potential energy of the beam
     * 
     * This method returns the strain energy stored in the beam due to deformation.
     * 
     * @param mparams Mechanical parameters
     * @param x Position vector
     * @return Potential energy
     */
    virtual double getPotentialEnergy(const sofa::core::MechanicalParams* mparams, 
                                    const DataVecCoord& x) const override;
    
    /**
     * @brief Returns the external radius of the beam
     * 
     * Utility method to access the beam's external radius.
     * 
     * @return The external radius of the beam
     */
    virtual Real getRadius();
};

#if !defined(SOFA_COMPONENT_FORCEFIELD_BASEBEAMHOOKELAWFORCEFIELD_CPP)
extern template class SOFA_COSSERAT_API BaseBeamHookeLawForceField<defaulttype::Vec3Types>;
extern template class SOFA_COSSERAT_API BaseBeamHookeLawForceField<defaulttype::Vec6Types>;
#endif

} // namespace sofa::component::forcefield

