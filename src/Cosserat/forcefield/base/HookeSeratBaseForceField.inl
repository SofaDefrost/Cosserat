#include <Cosserat/forcefield/base/HookeSeratBaseForceField.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/OptionsGroup.h>

namespace sofa::component::forcefield {

	template<typename DataTypes>
	HookeSeratBaseForceField<DataTypes>::HookeSeratBaseForceField() :
		ForceField<DataTypes>(), // Correct base class initialization
	d_crossSectionShape(initData(&d_crossSectionShape, {"circular", "rectangular", "hollow_circular"}, "crossSectionShape",
							 "shape of the cross-section. Can be: circular (solid circular), "
							 "rectangular (lengthY and lengthZ), or hollow_circular (tube with external "
							 "radius being radius and internal radius being innerRadius). Default is circular")),
		d_youngModulus(initData(&d_youngModulus, 1.0e9, "youngModulus",
								"Young Modulus describes the stiffness of the material")),
		d_poissonRatio(initData(&d_poissonRatio, 0.45, "poissonRatio",
								"Poisson Ratio describes the compressibility of the material")),
		d_length(initData(&d_length, "length", "The list of lengths of the different beam's sections.")),
		d_radius(initData(&d_radius, 1.0, "radius", "external radius of the cross section (if circular)")),
		d_innerRadius(
				initData(&d_innerRadius, 0.0, "innerRadius", "internal radius of the cross section (if circular)")),
		d_lengthY(initData(&d_lengthY, 1.0, "lengthY",
						   "side length of the cross section along local y axis (if rectangular)")),
		d_lengthZ(initData(&d_lengthZ, 1.0, "lengthZ",
						   "side length of the cross section along local z axis (if rectangular)")) {
		// Initialize computed properties
		m_crossSectionArea = 0.0;
		J = 0.0;
		Iy = 0.0;
		Iz = 0.0;
	}

	template<typename DataTypes>
	void HookeSeratBaseForceField<DataTypes>::init() {
		ForceField<DataTypes>::init(); // Correct base class call

		// Validate parameters before proceeding
		if (!validateParameters()) {
			msg_error("HookeSeratBaseForceField") << "Invalid parameters detected. Please check your configuration.";
			return;
		}

		reinit();
	}

	template<typename DataTypes>
	void HookeSeratBaseForceField<DataTypes>::reinit() {
		computeCrossSectionProperties();

		if (f_printLog.getValue()) {
			msg_info("HookeSeratBaseForceField") << "  ----------------------------------" ;
			printDebugInfo();
			msg_info("HookeSeratBaseForceField") << "  ----------------------------------" ;
		}
	}

	template<typename DataTypes>
	void HookeSeratBaseForceField<DataTypes>::computeCrossSectionProperties() {
		const std::string &shape = d_crossSectionShape.getValue().getSelectedItem();

		if (shape == "rectangular") {
			computeRectangularProperties();
		} else if (shape == "circular") {
			computeCircularProperties();
		} else if (shape == "hollow_circular") {
			computeHollowCircularProperties();
		} else {
			msg_error("HookeSeratBaseForceField") << "Unknown cross-section shape: " << shape;
		}
	}



	template<typename DataTypes>
void HookeSeratBaseForceField<DataTypes>::computeHollowCircularProperties() {
		const Real r = d_radius.getValue();
		const Real rInner = d_innerRadius.getValue();

		// Validate radii
		if (r <= 0) {
			msg_error("HookeSeratBaseForceField") << "External radius must be positive: " << r;
			return;
		}

		if (rInner <= 0) {
			msg_error("HookeSeratBaseForceField") << "Inner radius must be positive for hollow circular section: " << rInner;
			return;
		}

		if (rInner >= r) {
			msg_error("HookeSeratBaseForceField") << "Inner radius must be smaller than external radius: "
											<< "rInner=" << rInner << ", r=" << r;
			return;
		}

		const Real r2 = r * r;
		const Real r4 = r2 * r2;
		const Real rInner2 = rInner * rInner;
		const Real rInner4 = rInner2 * rInner2;

		// Cross-sectional area (annular area)
		m_crossSectionArea = M_PI * (r2 - rInner2);

		// Second moments of area (bending) - hollow circular section
		Iy = M_PI * (r4 - rInner4) / 4.0;
		Iz = Iy; // Circular symmetry

		// Torsional constant (polar moment of inertia) - hollow circular section
		J = M_PI * (r4 - rInner4) / 2.0;

		msg_info("HookeSeratBaseForceField") << "Hollow circular cross-section: r=" << r << ", rInner=" << rInner;
	}

	template<typename DataTypes>
void HookeSeratBaseForceField<DataTypes>::computeCircularProperties() {
		const Real r = d_radius.getValue();

		// Validate radius
		if (r <= 0) {
			msg_error("HookeSeratBaseForceField") << "External radius must be positive: " << r;
			return;
		}

		const Real r2 = r * r;
		const Real r4 = r2 * r2;

		// Cross-sectional area (solid circular)
		m_crossSectionArea = M_PI * r2;

		// Second moments of area (bending) - solid circular section
		Iy = M_PI * r4 / 4.0;
		Iz = Iy; // Circular symmetry

		// Torsional constant (polar moment of inertia) - solid circular section
		J = M_PI * r4 / 2.0;

		msg_info("HookeSeratBaseForceField") << "Solid circular cross-section: r=" << r;
	}

	template<typename DataTypes>
	void HookeSeratBaseForceField<DataTypes>::computeRectangularProperties() {
		const Real Ly = d_lengthY.getValue();
		const Real Lz = d_lengthZ.getValue();

		// Validate dimensions
		if (Ly <= 0 || Lz <= 0) {
			msg_error("HookeSeratBaseForceField") << "Rectangular dimensions must be positive: "
											<< "Ly=" << Ly << ", Lz=" << Lz;
			return;
		}

		// Cross-sectional area
		m_crossSectionArea = Ly * Lz;

		// Second moments of area
		Iy = (Ly * Lz * Lz * Lz) / 12.0; // About y-axis
		Iz = (Lz * Ly * Ly * Ly) / 12.0; // About z-axis

		// Torsional constant for rectangular section
		// Using Saint-Venant's theory approximation
		const Real a = std::max(Ly, Lz); // Larger dimension
		const Real b = std::min(Ly, Lz); // Smaller dimension
		const Real ratio = b / a;

		// Approximation formula for torsional constant
		Real beta;
		if (ratio >= 1.0) {
			beta = 1.0 / 3.0; // Square section
		} else {
			// Approximation for rectangular sections
			beta = (1.0 / 3.0) * (1.0 - 0.63 * ratio + 0.052 * ratio * ratio * ratio * ratio * ratio);
		}

		J = beta * a * b * b * b; // CORRECTED: proper torsional constant

		msg_info("HookeSeratBaseForceField") << "Rectangular cross-section: Ly=" << Ly << ", Lz=" << Lz;
	}

	template<typename DataTypes>
	bool HookeSeratBaseForceField<DataTypes>::validateParameters() const {
		// Check material properties
		if (d_youngModulus.getValue() <= 0) {
			msg_error("HookeSeratBaseForceField") << "Young's modulus must be positive: " << d_youngModulus.getValue();
			return false;
		}

		const Real nu = d_poissonRatio.getValue();
		if (nu <= -1.0 || nu >= 0.5) {
			msg_error("HookeSeratBaseForceField") << "Poisson's ratio must be in range (-1, 0.5): " << nu;
			return false;
		}

		// Check geometric parameters
		const std::string &shape = d_crossSectionShape.getValue().getSelectedItem();

		if (shape == "circular") {
			if (d_radius.getValue() <= 0) {
				msg_error("HookeSeratBaseForceField") << "Radius must be positive: " << d_radius.getValue();
				return false;
			}
		} else if (shape == "hollow_circular") {
			if (d_radius.getValue() <= 0) {
				msg_error("HookeSeratBaseForceField") << "External radius must be positive: " << d_radius.getValue();
				return false;
			}
			if (d_innerRadius.getValue() <= 0) {
				msg_error("HookeSeratBaseForceField") << "Inner radius must be positive for hollow circular section: " << d_innerRadius.getValue();
				return false;
			}
			if (d_innerRadius.getValue() >= d_radius.getValue()) {
				msg_error("HookeSeratBaseForceField") << "Inner radius must be smaller than external radius: "
												<< "rInner=" << d_innerRadius.getValue() << ", r=" << d_radius.getValue();
				return false;
			}
		} else if (shape == "rectangular") {
			if (d_lengthY.getValue() <= 0 || d_lengthZ.getValue() <= 0) {
				msg_error("HookeSeratBaseForceField") << "Rectangular dimensions must be positive: "
												<< "Ly=" << d_lengthY.getValue() << ", Lz=" << d_lengthZ.getValue();
				return false;
			}
		}

		return true;
	}

	template<typename DataTypes>
	bool HookeSeratBaseForceField<DataTypes>::isValidConfiguration() const {
		return validateParameters() && (m_crossSectionArea > 0) && (J > 0) && (Iy > 0) && (Iz > 0);
	}

	template<typename DataTypes>
	void HookeSeratBaseForceField<DataTypes>::printDebugInfo() const {
		msg_info("HookeSeratBaseForceField") << "  ----------------------------------" ;
		msg_info("HookeSeratBaseForceField") << "Cross-section properties:";
		msg_info("HookeSeratBaseForceField") << "  Shape: " << d_crossSectionShape.getValue().getSelectedItem();
		msg_info("HookeSeratBaseForceField") << "  Area: " << m_crossSectionArea;
		msg_info("HookeSeratBaseForceField") << "  J (torsion): " << J;
		msg_info("HookeSeratBaseForceField") << "  Iy (bending): " << Iy;
		msg_info("HookeSeratBaseForceField") << "  Iz (bending): " << Iz;
		msg_info("HookeSeratBaseForceField") << "  E (Young): " << d_youngModulus.getValue();
		msg_info("HookeSeratBaseForceField") << "  ν (Poisson): " << d_poissonRatio.getValue();
		msg_info("HookeSeratBaseForceField") << "  G (Shear): " << getShearModulus();
		msg_info("HookeSeratBaseForceField") << "  ----------------------------------" ;
	}

} // namespace sofa::component::forcefield
