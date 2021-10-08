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

#include "BeamPlasticLawForceField.h"

using sofa::helper::ReadAccessor;
using sofa::helper::WriteAccessor;
using sofa::core::VecCoordId;

namespace sofa::component::forcefield
{

template<typename DataTypes>
BeamPlasticLawForceField<DataTypes>::BeamPlasticLawForceField() :
    Inherit1(),
    d_initialYieldStress(initData(&d_initialYieldStress, "initialYieldStress",
        "Yield stress of the considered material, prior to any elastic deformation")),
    d_plasticModulus(initData(&d_plasticModulus, "plasticModulus",
        "Approximation of the plastic modulus as a constant. Can be deduced from a generic law such as Ramberg-Osgood's"))
{
    if (d_initialYieldStress.getValue() < 0)
    {
        msg_error() << "yield Stress should be positive. Please provide a positive yield"
                    << "stress value for the considered material";
    }
}


template<typename DataTypes>
BeamPlasticLawForceField<DataTypes>::~BeamPlasticLawForceField()
{
}


template<typename DataTypes>
void BeamPlasticLawForceField<DataTypes>::reinit()
{
    Inherit1::reinit();

    // Initialisation of the lastStrain field with the rest strain
    m_lastStrain = this->mstate->read(VecCoordId::restPosition())->getValue();

    size_t nbSections = this->d_length.getValue().size();

    // Initialisation of the prevStress field with 0
    m_prevStress.clear();
    m_prevStress.resize(nbSections, Vec3());

    // Initialisaiton of the section mechanical states to ELASTIC
    m_sectionMechanicalStates.clear();
    m_sectionMechanicalStates.resize(nbSections, MechanicalState::ELASTIC);

    // Initialisaiton of the tangent stiffness matrices with the elastic stiffness matrices
    m_Kt_sectionList.clear();
    if (!this->d_varianteSections.getValue())
    {
        m_Kt_sectionList.resize(nbSections, this->m_K_section);
    }
    else
    {
        for (unsigned int i = 0; i < nbSections; i++)
            m_Kt_sectionList.push_back(this->m_K_sectionList[i]);
    }

    // Computation of the generalised Hooke's law
    // As we are working with only 3 components of the strain tensor,
    // the generalised Hooke's law is reduced to a 3x3 diagonal matrix
    // (instead of a 4th order tensor represented as a 9x9 matrix)
    if (!this->d_varianteSections.getValue()) {

        Real E = this->d_youngModulus.getValue();
        Real nu = this->d_poissonRatio.getValue();

        m_genHookesLaw[0][0] = m_genHookesLaw[1][1] = m_genHookesLaw[2][2] = E / (1 + nu);
    }
    else
    {
        m_genHookesLawList.clear();
        for (unsigned int i = 0; i < nbSections; i++)
        {
            Mat33 _genHookesLaw;
            Real E = this->d_youngModulusList.getValue()[i];
            Real nu = this->d_poissonRatioList.getValue()[i];

            _genHookesLaw[0][0] = _genHookesLaw[1][1] = _genHookesLaw[2][2] = E / (1 + nu);
            m_genHookesLawList.push_back(_genHookesLaw);
        }
    }

    // Initialisation of plasticity parameters
    m_backStress.clear();
    m_backStress.resize(nbSections, Vec3());

    m_yieldStress.clear();
    m_yieldStress.resize(nbSections, d_initialYieldStress.getValue());

    //By default, no plastic deformation => no history
    m_plasticStrain.clear();
    m_plasticStrain.resize(nbSections, Vec3());

    m_effectivePlasticStrain.clear();
    m_effectivePlasticStrain.resize(nbSections, 0);

    // Initialisation of the comparison threshold for stress tensor norms to 0.
    // Plasticity computation requires to basically compare stress tensor norms to 0.
    // As stress norm values can vary of several orders of magnitude, depending on the
    // considered materials and/or applied forces, this comparison has to be carried out
    // carefully.
    // The idea here is to use the initialYieldStress of the material, and the
    // available precision limit (e.g. std::numeric_limits<double>::epsilon()).
    // We rely on the value of the initial Yield stress, as we can expect plastic
    // deformation to occur inside a relatively small intervl of stresses around this value.
    const int orderOfMagnitude = d_initialYieldStress.getValue(); //Should use std::abs, but d_initialYieldStress > 0
    m_stressComparisonThreshold = std::numeric_limits<double>::epsilon() * orderOfMagnitude;
}


template< class DataTypes>
typename BeamPlasticLawForceField<DataTypes>::Vec3 BeamPlasticLawForceField<DataTypes>::vonMisesGradient(const Vec3& stressTensor)
{
    // NB: this gradient represent the normal to the yield surface
    // in case the Von Mises yield criterion is used. The norm of the
    // gradient is sqrt(3/2): it has to be multiplied by sqrt(2/3)
    // to give the unit normal to the yield surface

    // Deviatoric-based computation
    Real sigmaEq = equivalentStress(stressTensor);
    // In the following computation, we should rigorously use the deviatoric
    // counterpart of stressTensor, instead of stressTensor directly. Here,
    // the Cosserat model we use reduces only to the non-diagonal components
    // of the stress tensor. Therefore the reduced stress-tensor considered
    // here is deviatoric.
    if (sigmaEq < m_stressComparisonThreshold) // We consider that sigmaEq = 0
        return Vec3();
    else
        return (3 / (2 * sigmaEq)) * stressTensor;
}


template<typename DataTypes>
void BeamPlasticLawForceField<DataTypes>::computeStressIncrement(unsigned int sectionId,
                                                                 const Coord& strainIncrement,
                                                                 Vec3& newStressPoint)
{
    /// This method implements the radial return algorithm, as in "Numerical Implementation of
    /// Constitutive models: Rate-independent Deviatoric Plasticity", T.J.R. Hugues, 1984.
    /// The idea is to compute the stress increment in two steps : a purely elastic step, in
    /// which all deformation is considered elastic, and a plastic correction step, is
    /// deformation was actually important enough to generate plasticity.
    /// The plasticity model used in the computation is a Von Mises-Hill plasticity with
    /// linear mixed hardening.
    /// NB: we consider that the yield function and the plastic flow are equal (f=g). This
    /// corresponds to an associative flow rule (for plasticity).

    const Mat33& C = m_genHookesLaw;

    //// First step = computation of the elastic predictor, as if deformation was entirely elastic
    const MechanicalState mechanicalState = m_sectionMechanicalStates[sectionId];

    Vec3 elasticIncrement = C * strainIncrement;
    Vec3 elasticPredictorStress = m_prevStress[sectionId] + elasticIncrement;

    const Vec3& backStress = m_backStress[sectionId];
    const Real yieldStress = m_yieldStress[sectionId];

    if (vonMisesYield(elasticPredictorStress, backStress, yieldStress) < m_stressComparisonThreshold)
    {
        // The segment is in elastic state: the back stress and yield stress
        // remain constant, and the new stress is equal to the trial stress.
        newStressPoint = elasticPredictorStress;

        // If the segment was initially plastic, we update its mechanical state
        if (mechanicalState == MechanicalState::PLASTIC)
            m_sectionMechanicalStates[sectionId] = MechanicalState::POSTPLASTIC;
    }
    else
    {
        // If the segment was initially elastic, we update its mechanical state
        if (mechanicalState == MechanicalState::POSTPLASTIC || mechanicalState == MechanicalState::ELASTIC)
            m_sectionMechanicalStates[sectionId] = MechanicalState::PLASTIC;

        // /!\ We here consider that we obtain a deviatoric stress tensor, as the diagonal components of a
        // complete stress tensor, representing the axial stresses, are ignored in the model.
        Vec3 shiftedDeviatoricElasticPredictor = elasticPredictorStress - backStress;

        // Gradient of the Von Mises yield function is colinear to the deviatoric stress tensor.
        // Thus we can compute the yield surface normal using the deviatoric stress.
        // For the Von Mises yield function, the normal direction to the yield surface doesn't
        // change between the elastic predictor and it's projection on the yield surface
        Real shiftDevElasticPredictorNorm = correctedNorm(shiftedDeviatoricElasticPredictor);
        Vec3 N = shiftedDeviatoricElasticPredictor / shiftDevElasticPredictorNorm;

        // Indicates the proportion of Kinematic vs isotropic hardening. beta=0 <=> kinematic, beta=1 <=> isotropic
        const Real beta = 0.5;

        Real E, nu = 0;
        if (!this->d_varianteSections.getValue())
        {
            E = this->d_youngModulus.getValue();
            nu = this->d_poissonRatio.getValue();
        }
        else
        {
            E = this->d_youngModulusList.getValue()[sectionId];
            nu = this->d_poissonRatioList.getValue()[sectionId];
        }
        const Real mu = E / (2 * (1 + nu)); // Lame coefficient

        // Plastic modulus
        const Real H = d_plasticModulus.getValue();

        // Computation of the plastic multiplier
        const double sqrt2 = helper::rsqrt(2.0);
        const double sqrt3 = helper::rsqrt(3.0);
        const double sqrt6 = sqrt2 * sqrt3;
        Real plasticMultiplier = (shiftDevElasticPredictorNorm - (sqrt2 / sqrt3) * yieldStress) / (mu * sqrt6 * (1 + H / (3 * mu)));

        newStressPoint = elasticPredictorStress - sqrt6 * mu * plasticMultiplier * N;

        // Updating plastic variables
        Real newYieldStress = yieldStress + beta * H * plasticMultiplier;
        m_yieldStress[sectionId] = newYieldStress;
        Vec3 newBackStress = backStress + (sqrt2 / sqrt3) * (1 - beta) * H * plasticMultiplier * N;
        m_backStress[sectionId] = newBackStress;

        Vec3 plasticStrainIncrement = (sqrt3 / sqrt2) * plasticMultiplier * N;
        m_plasticStrain[sectionId] += plasticStrainIncrement;
        m_effectivePlasticStrain[sectionId] += plasticMultiplier;
    }

    // Updating the stress value, wether elastic or plastic
    m_prevStress[sectionId] = newStressPoint;
}


template<typename DataTypes>
void BeamPlasticLawForceField<DataTypes>::updateTangentStiffness(unsigned int sectionId)
{
    Mat33 Kt = Mat33();

    const Mat33& C = m_genHookesLaw;
    Real E, nu = 0;
    if (!this->d_varianteSections.getValue())
    {
        E = this->d_youngModulus.getValue();
        nu = this->d_poissonRatio.getValue();
    }
    else
    {
        E = this->d_youngModulusList.getValue()[sectionId];
        nu = this->d_poissonRatioList.getValue()[sectionId];
    }

    Vec3 currentStressPoint = m_prevStress[sectionId]; //Updated in computeStressIncrement

    Real H = d_plasticModulus.getValue();

    Mat33 Cep = Mat33();
    // Cep
    Vec3 gradient = vonMisesGradient(currentStressPoint);

    if (correctedNorm(gradient) < m_stressComparisonThreshold
        || m_sectionMechanicalStates[sectionId] != MechanicalState::PLASTIC)
        Cep = C; //TO DO: is that correct ?
    else
    {
        Vec3 N = helper::rsqrt(2.0 / 3.0) * gradient; // Normal to the yield surface
        Vec3 CN = C * N;
        //Conversion to matrix, TO DO: better way ? Eigen ? Use only matrices ?
        Mat<3, 1, Real> matCN = Mat<3, 1, Real>();
        for (unsigned int i = 0; i < CN.size(); i++)
            matCN[i][0] = CN[i];
        // NtC = (NC)t because of C symmetry
        Cep = C - (2 * matCN * matCN.transposed()) / (2.0 * N * CN + (2.0 / 3.0) * H); //TO DO: check that * operator is actually dot product

        // /!\ Warning on the computation above /!\
        // Terms CNNtC in the numerator and NtCN in the denominator are multiplied by 2
        // in order to account for the fact that we are using a reduced notation.
        // In the same way that the generalised Hooke's law we use in m_genHookesLaw
        // differs from the 'complete' 9x9 generalised Hooke's Law components by a
        // factor 2, we multiply Cep components by 2 to account for the fact that
        // multiplication of this matrix by a vector representing a symmetric tensor
        // makes each non-diagonal terms appear twice.
        // NB: this takes into account the fact that C is already expressed here with
        // a factor 2.
        // /!\ This should be used with serious caution in all computations, as it
        // always has to be coherent with the actual complete tensor computation.
        // TO DO: better way to implement this ?
    }

    // Integration step, consisting only in multiplication by the element volume, as
    // Cep is considered constant over a segment.
    Real L = this->d_length.getValue()[sectionId];
    Real A = this->m_crossSectionArea;
    Real V = L * A; // Assuming the initial volume is conserved during deformation
    Kt = V * Cep;

    // Update the beam tangent stiffness
    m_Kt_sectionList[sectionId] = Kt;
}


template<typename DataTypes>
void BeamPlasticLawForceField<DataTypes>::addForce(const MechanicalParams* mparams,
                                                   DataVecDeriv& d_f,
                                                   const DataVecCoord& d_x,
                                                   const DataVecDeriv& d_v)
{
    SOFA_UNUSED(d_v);
    SOFA_UNUSED(mparams);

    if (!this->getMState()) {
        msg_info("BeamPlasticLawForceField") << "No Mechanical State found, no force will be computed..." << "\n";
        this->compute_df = false;
        return;
    }

    VecDeriv& f = *d_f.beginEdit();
    const VecCoord& x = d_x.getValue();

    f.resize(x.size());

    if (x.size() != this->d_length.getValue().size()) {
        msg_warning("BeamPlasticLawForceField") << " length should have the same size as x..." << "\n";
        this->compute_df = false;
        return;
    }

    for (unsigned int i = 0; i < x.size(); i++)
    {
        // Strain is considered constant over each section, so Gaussian integration is not required
        const Coord strainIncrement = x[i] - m_lastStrain[i];

        // Computation of the stress using the radial return algorithm
        Vec3 newStressPoint = Vec3();
        computeStressIncrement(i, strainIncrement, newStressPoint);

        // If the beam is in plastic state, we update the tangent stiffness matrix
        updateTangentStiffness(i);

        // Computation of internal forces from stress
        // As stress and strain are uniform over the segment, spatial integration reduces to the segment volume
        Real L = this->d_length.getValue()[i];
        Real A = this->m_crossSectionArea;
        Real V = L * A; // Assuming the initial volume is conserved during deformation
        f[i] = V * m_prevStress[i]; // Previous stress has been updated in computeStressIncrement
    }

    // Save the current strain as a record for the next time step.
    //TO DO: check if this is copy operator
    m_lastStrain = x;

    d_f.endEdit();
}


template<typename DataTypes>
void BeamPlasticLawForceField<DataTypes>::addDForce(const MechanicalParams* mparams,
                                                    DataVecDeriv& d_df,
                                                    const DataVecDeriv& d_dx)
{
    if (!this->compute_df)
        return;

    WriteAccessor< DataVecDeriv > df = d_df;
    ReadAccessor< DataVecDeriv > dx = d_dx;
    Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    df.resize(dx.size());

    for (unsigned int i = 0; i < dx.size(); i++)
        df[i] -= (m_Kt_sectionList[i] * dx[i]) * kFactor * this->d_length.getValue()[i];
}


template<typename DataTypes>
void BeamPlasticLawForceField<DataTypes>::addKToMatrix(const MechanicalParams* mparams,
                                                       const MultiMatrixAccessor* matrix)
{
    MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
    BaseMatrix* mat = mref.matrix;
    unsigned int offset = mref.offset;
    Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());

    const VecCoord& pos = this->mstate->read(core::ConstVecCoordId::position())->getValue();

    for (unsigned int n = 0; n < pos.size(); n++)
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                mat->add(offset + i + 3 * n, offset + j + 3 * n, -kFact * m_Kt_sectionList[n][i][j] * this->d_length.getValue()[n]);
}


template< class DataTypes>
typename BeamPlasticLawForceField<DataTypes>::Real BeamPlasticLawForceField<DataTypes>::equivalentStress(const Vec3& stressTensor)
{
    // Direct computation of the equivalent stress. /!\ Only 3 components of
    // the strain tensor are considered in this Cosserat model (representing
    // torsion and bending). Consequently, we consider only the 3 components
    // of the stress tensor that would result from these 3 strain components
    // in a purely elastic deformation (in which the stress tensor would be
    // obtained from strain using Hooke's law).
    // This is a priori an approximation as the relation between stress and
    // strain during plastic deformation is more complexe.

    Real sigmaYZ = stressTensor[0];
    Real sigmaXZ = stressTensor[1];
    Real sigmaXY = stressTensor[2];

    return helper::rsqrt( 3.0 * (sigmaYZ * sigmaYZ + sigmaXZ * sigmaXZ + sigmaXY * sigmaXY) );
}


template< class DataTypes>
typename BeamPlasticLawForceField<DataTypes>::Real BeamPlasticLawForceField<DataTypes>::vonMisesYield(const Vec3& stressTensor,
                                                                                                      const Vec3& backStress,
                                                                                                      const Real yieldStress)
{
    return equivalentStress(stressTensor - backStress) - yieldStress;
}


template< class DataTypes>
typename BeamPlasticLawForceField<DataTypes>::Real BeamPlasticLawForceField<DataTypes>::correctedNorm(const Vec3& tensor2)
{
    Real tensorYZ = tensor2[0];
    Real tensorXZ = tensor2[1];
    Real tensorXY = tensor2[2];

    return helper::rsqrt(2*tensorYZ*tensorYZ + 2*tensorXZ*tensorXZ + 2*tensorXY*tensorXY);
}


template< class DataTypes>
const sofa::type::vector<typename BeamPlasticLawForceField<DataTypes>::MechanicalState>& BeamPlasticLawForceField<DataTypes>::getSectionMechanicalStates()
{
    return m_sectionMechanicalStates;
}

} // sofa::component::forcefield
