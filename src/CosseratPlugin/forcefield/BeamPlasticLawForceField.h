/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
#pragma once
#include <CosseratPlugin/config.h>

#include <CosseratPlugin/forcefield/BeamHookeLawForceField.h>

namespace sofa::component::forcefield
{

using sofa::type::vector;
using sofa::type::Vec;

template<typename DataTypes>
class BeamPlasticLawForceField : public BeamHookeLawForceField<DataTypes>
{

public:

    SOFA_CLASS(SOFA_TEMPLATE(BeamPlasticLawForceField, DataTypes),
        SOFA_TEMPLATE(BeamHookeLawForceField, DataTypes));

    typedef typename DataTypes::Real     Real;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord;
    typedef typename BeamHookeLawForceField<DataTypes>::DeformationRegime DeformationRegime;

    typedef Data<VecCoord>    DataVecCoord;
    typedef Data<VecDeriv>    DataVecDeriv;

    typedef Vec<3, Real>                Vec3;
    typedef Mat<3, 3, Real>             Mat33;


public:

    BeamPlasticLawForceField();
    virtual ~BeamPlasticLawForceField();

    ////////////////////////// Inherited from BaseObject /////////////////////////
    void reinit();
    ///////////////////////////////////////////////////////////////////////////

    ////////////////////////// Inherited from ForceField /////////////////////////
    void addForce(const MechanicalParams* mparams,
                  DataVecDeriv& f,
                  const DataVecCoord& x,
                  const DataVecDeriv& v) override;

    void addDForce(const MechanicalParams* mparams,
                   DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    void addKToMatrix(const MechanicalParams* mparams,
                      const MultiMatrixAccessor* matrix) override;
    ///////////////////////////////////////////////////////////////////////////

    // Inherited from BeamHookeLawForceField
    bool isPlastic() const;

protected:

    Data<vector<Real>> d_initialYieldStresses;

    // TO DO: Plastic modulus should not be a constant. Is this approximation relevant ?
    // Otherwise a generic law such as Ramberg-Osgood should be use. (Cf Plugin BeamPlastic)
    Data<vector<Real>> d_plasticModuli; ///Linearisation coefficient for the plastic behaviour law

    /// Coefficients to determine the porportion of isotropic and kinematic hardening
    /// 0 = purely kinematic
    /// 1 = purely isotropic
    Data<vector<Real>> d_mixedHardeningCoefficients;


    // Threshold used to compare stress tensor norms to 0. See detailed explanation
    // at the computation of the threshold in the init() method.
    Real m_stressComparisonThreshold;

    /// List of the strain coordinates associated to the length sections, at the previous time step
    VecCoord m_lastStrain;
    /// List of the stress associated to each length sections, at the previous time step
    vector<Vec3> m_prevStress;
    /// List of tangent stiffness matrices, taking into account the plastic deformation of the length sections
    vector<Mat33> m_Kt_sectionList;

    //Plasticity history variables
    vector<Vec3> m_plasticStrain;
    vector<Real> m_effectivePlasticStrain;

    vector<Vec3> m_backStress; /// Centre of the yield surface, in stress space (one for each length section)
    // NB: the distinction between d_initialYieldStresses and m_yieldStress allows to handle the evolution of
    // the yield stresses internally, without accessing any user data. The evolution of yield stress should
    // depend entirely on the plasticity method, and it doesn't make sense for a user to have access to it
    // (left aside for the parameter initialisation).
    vector<Real> m_yieldStress; /// Elastic limit, varying if plastic deformation occurs (one for each length section)

    /// Generalised Hooke's law, reduced to the considered strain and stress components
    Mat33 m_genHookesLaw;
    /// List of generalised Hooke's law if each segment has a specific Young's modulus / Poisson ratio
    vector<Mat33> m_genHookesLawList;

    //----- Implementation of the Von Mises yield function -----//

    /// Computes the equivalent stress from a tensor
    Real equivalentStress(const Vec3& stressTensor);
    /// Evaluates the Von Mises yield function for given stress tensor and yield stress
    Real vonMisesYield(const Vec3& stressTensor, const Vec3& backStress, const Real yieldStress);
    /// Computes the Von Mises yield function gradient, for a given stress tensor
    Vec3 vonMisesGradient(const Vec3& stressTensor);

    //----- Correction methods for reduced notation -----//

    /// Computes the norm of a Vec3 reduced representation of a 2nd-order tensor
    Real correctedNorm(const Vec3& tensor2);

    /// Computes stress increment on a single point of space, from previous stress and current strain
    void computeStressIncrement(unsigned int sectionId, const Coord& strainIncrement,
                                Vec3& newStressPoint);

    void updateTangentStiffness(unsigned int sectionId);

};

}  // sofa::component::forcefield
