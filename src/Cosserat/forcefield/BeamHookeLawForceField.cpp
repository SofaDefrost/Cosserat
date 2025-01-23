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
#define SOFA_COSSERAT_CPP_BeamHookeLawForceField
#include <Cosserat/forcefield/BeamHookeLawForceField.inl>

#include <sofa/core/ObjectFactory.h>

namespace sofa::component::forcefield {

template <> void BeamHookeLawForceField<defaulttype::Vec6Types>::reinit() {
  // Precompute and store values
  Real Iy, Iz, J, A;
  if (d_crossSectionShape.getValue().getSelectedItem() ==
      "rectangular") // rectangular cross-section
  {
    Real Ly = d_lengthY.getValue();
    Real Lz = d_lengthZ.getValue();

    const Real LyLzLzLz = Ly * Lz * Lz * Lz;
    const Real LzLyLyLy = Lz * Ly * Ly * Ly;

    Iy = LyLzLzLz / 12.0;
    Iz = LzLyLyLy / 12.0;
    J = Iy + Iz;
    A = Ly * Lz;

  } else // circular cross-section
  {
    msg_info() << "Cross section shape."
               << d_crossSectionShape.getValue().getSelectedItem();

    Real r = d_radius.getValue();
    Real rInner = d_innerRadius.getValue();
    const Real r4 = r * r * r * r;
    const Real rInner4 = rInner * rInner * rInner * rInner;

    Iy = M_PI * (r4 - rInner4) / 4.0;
    Iz = Iy;
    J = Iy + Iz;
    A = M_PI * (r * r - rInner4);
  }
  m_crossSectionArea = A;

  if (d_useInertiaParams.getValue()) {
    msg_info("BeamHookeLawForceField")
        << "Pre-calculated inertia parameters are used for the computation of the stiffness matrix.";
    m_K_section66[0][0] = d_GI.getValue();
    m_K_section66[1][1] = d_EIy.getValue();
    m_K_section66[2][2] = d_EIz.getValue();
    m_K_section66[3][3] = d_EA.getValue();
    m_K_section66[4][4] = d_GA.getValue();
    m_K_section66[5][5] = d_GA.getValue();
  } else {
    // Pour du vec 6, on a  _m_K =diag([G*J E*Iy E*Iz E*A G*As G*As]); % stifness matrix
    Real E = d_youngModulus.getValue();
    Real G = E / (2.0 * (1.0 + d_poissonRatio.getValue()));
    // Translational Stiffness (X, Y, Z directions):
    //  Gs * J(i): Shearing modulus times the second moment of the area (torsional stiffness). E * Js(i): Young's modulus times the second moment of the area (bending stiffness).
    m_K_section66[0][0] = G * J;
    m_K_section66[1][1] = E * Iy;
    m_K_section66[2][2] = E * Iz;
    // Rotational Stiffness (X, Y, Z directions):
    // E * A: Young's modulus times the cross-sectional area (axial stiffness).
    // Gs * A: Shearing modulus times the cross-sectional area (shearing stiffness).
    m_K_section66[3][3] = E * A;
    m_K_section66[4][4] = G * A;
    m_K_section66[5][5] = G * A;
  }
}

template <>
void BeamHookeLawForceField<defaulttype::Vec6Types>::addForce(
    const MechanicalParams *mparams, DataVecDeriv &d_f, const DataVecCoord &d_x,
    const DataVecDeriv &d_v) {
  SOFA_UNUSED(d_v);
  SOFA_UNUSED(mparams);

  if (!this->getMState()) {
    msg_info("BeamHookeLawForceField")
        << "No Mechanical State found, no force will be computed..." << "\n";
    compute_df = false;
    return;
  }
  VecDeriv &f = *d_f.beginEdit();
  const VecCoord &x = d_x.getValue();
  // get the rest position (for non straight shape)
  const VecCoord &x0 =
      this->mstate->read(VecCoordId::restPosition())->getValue();

  f.resize(x.size());
  unsigned int sz = d_length.getValue().size();
  if (x.size() != sz) {
    msg_warning("BeamHookeLawForceField")
        << " length : " << sz << "should have the same size as x... "
        << x.size() << "\n";
    compute_df = false;
    return;
  }
  for (unsigned int i = 0; i < x.size(); i++)
    f[i] -= (m_K_section66 * (x[i] - x0[i])) * d_length.getValue()[i];

  d_f.endEdit();
}

template <>
void BeamHookeLawForceField<defaulttype::Vec6Types>::addDForce(
    const MechanicalParams *mparams, DataVecDeriv &d_df,
    const DataVecDeriv &d_dx) {
  if (!compute_df)
    return;

  WriteAccessor<DataVecDeriv> df = d_df;
  ReadAccessor<DataVecDeriv> dx = d_dx;
  Real kFactor = (Real)mparams->kFactorIncludingRayleighDamping(
      this->rayleighStiffness.getValue());

  df.resize(dx.size());
  for (unsigned int i = 0; i < dx.size(); i++)
    df[i] -= (m_K_section66 * dx[i]) * kFactor * d_length.getValue()[i];
}

template <>
void BeamHookeLawForceField<defaulttype::Vec6Types>::addKToMatrix(
    const MechanicalParams *mparams, const MultiMatrixAccessor *matrix) {
  MultiMatrixAccessor::MatrixRef mref = matrix->getMatrix(this->mstate);
  BaseMatrix *mat = mref.matrix;
  unsigned int offset = mref.offset;
  Real kFact = (Real)mparams->kFactorIncludingRayleighDamping(
      this->rayleighStiffness.getValue());

  const VecCoord &pos =
      this->mstate->read(core::ConstVecCoordId::position())->getValue();
  for (unsigned int n = 0; n < pos.size(); n++) {
    for (unsigned int i = 0; i < 6; i++)
      for (unsigned int j = 0; j < 6; j++)
        mat->add(offset + i + 6 * n, offset + j + 6 * n,
                 -kFact * m_K_section66[i][j] * d_length.getValue()[n]);
  }
}

using namespace sofa::defaulttype;

template class BeamHookeLawForceField<Vec3Types>;
template class BeamHookeLawForceField<Vec6Types>;

} // namespace Cosserat::forcefield

// namespace sofa
namespace Cosserat
{

void registerBeamHookeLawForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(
        sofa::core::ObjectRegistrationData(
            "This component is used to compute internal stress for torsion (along x) and bending (y and z)")
            .add<sofa::component::forcefield::BeamHookeLawForceField<sofa::defaulttype::Vec3Types>>(true)
            .add<sofa::component::forcefield::BeamHookeLawForceField<sofa::defaulttype::Vec6Types>>());
}

}
