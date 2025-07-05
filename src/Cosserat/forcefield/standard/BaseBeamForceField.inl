#include <Cosserat/forcefield/standard/BaseBeamForceField.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/OptionsGroup.h> // ??
#include <sofa/linearalgebra/FullVector.h>

namespace sofa::component::forcefield {

	template<typename DataTypes>
	BaseBeamForceField<DataTypes>::BaseBeamForceField() :
		d_crossSectionShape(initData(&d_crossSectionShape, {"circular", "rectangular"}, "crossSectionShape",
									 "shape of the cross-section. Can be: circular (tube with external "
									 "radius being radius and internal radius being innerRadius ) or "
									 "rectangular (lengthY and lengthZ) . Default is circular")),
		d_youngModulus(initData(&d_youngModulus, 1.0e9, "youngModulus",
								"Young Modulus describes the stiffness of the material")),
		d_poissonRatio(initData(&d_poissonRatio, 0.45, "poissonRatio",
								"poisson Ratio describes the compressibility of the material")),
		d_length(initData(&d_length, "length", "The list of lengths of the different beam's sections.")),
		d_radius(initData(&d_radius, 1.0, "radius", "external radius of the cross section (if circular)")),
		d_innerRadius(
				initData(&d_innerRadius, 0.0, "innerRadius", "internal radius of the cross section (if circular)")),
		d_lengthY(initData(&d_lengthY, 1.0, "lengthY",
						   "side length of the cross section along local y axis "
						   "(if rectangular)")),
		d_lengthZ(initData(&d_lengthZ, 1.0, "lengthZ",
						   "side length of the cross section along local z axis "
						   "(if rectangular)")) {}

template<typename DataTypes>
void BaseBeamForceField<DataTypes>::init() {
	Inherit1::init();
		reinit();
	}

	template<typename DataTypes>
	void BaseBeamForceField<DataTypes>::reinit() {
		// Precompute and store inertia values
		Real A;
		if (d_crossSectionShape.getValue().getSelectedItem() == "rectangular") // rectangular cross-section
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
			msg_info() << "Cross section shape." << d_crossSectionShape.getValue().getSelectedItem();

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
	}

} // namespace sofa::component::forcefield
