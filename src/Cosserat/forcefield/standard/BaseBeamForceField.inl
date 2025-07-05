#include <Cosserat/forcefield/standard/BaseBeamForceField.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/OptionsGroup.h>

namespace sofa::component::forcefield {

	template<typename DataTypes>
	BaseBeamForceField<DataTypes>::BaseBeamForceField() :
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
	void BaseBeamForceField<DataTypes>::init() {
		ForceField<DataTypes>::init(); // Correct base class call

		// Validate parameters before proceeding
		if (!validateParameters()) {
			msg_error("BaseBeamForceField") << "Invalid parameters detected. Please check your configuration.";
			return;
		}

		reinit();
	}

	template<typename DataTypes>
	void BaseBeamForceField<DataTypes>::reinit() {
		computeCrossSectionProperties();

		if (f_printLog.getValue()) {
			printDebugInfo();
		}
	}

	template<typename DataTypes>
	void BaseBeamForceField<DataTypes>::computeCrossSectionProperties() {
		const std::string &shape = d_crossSectionShape.getValue().getSelectedItem();

		if (shape == "rectangular") {
			computeRectangularProperties();
		} else if (shape == "circular") {
			computeCircularProperties();
		} else if (shape == "hollow_circular") {
			computeHollowCircularProperties();
		} else {
			msg_error("BaseBeamForceField") << "Unknown cross-section shape: " << shape;
		}
	}



	template<typename DataTypes>
void BaseBeamForceField<DataTypes>::computeHollowCircularProperties() {
		const Real r = d_radius.getValue();
		const Real rInner = d_innerRadius.getValue();

		// Validate radii
		if (r <= 0) {
			msg_error("BaseBeamForceField") << "External radius must be positive: " << r;
			return;
		}

		if (rInner <= 0) {
			msg_error("BaseBeamForceField") << "Inner radius must be positive for hollow circular section: " << rInner;
			return;
		}

		if (rInner >= r) {
			msg_error("BaseBeamForceField") << "Inner radius must be smaller than external radius: "
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

		msg_info("BaseBeamForceField") << "Hollow circular cross-section: r=" << r << ", rInner=" << rInner;
	}

	template<typename DataTypes>
void BaseBeamForceField<DataTypes>::computeCircularProperties() {
		const Real r = d_radius.getValue();

		// Validate radius
		if (r <= 0) {
			msg_error("BaseBeamForceField") << "External radius must be positive: " << r;
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

		msg_info("BaseBeamForceField") << "Solid circular cross-section: r=" << r;
	}

	template<typename DataTypes>
	void BaseBeamForceField<DataTypes>::computeRectangularProperties() {
		const Real Ly = d_lengthY.getValue();
		const Real Lz = d_lengthZ.getValue();

		// Validate dimensions
		if (Ly <= 0 || Lz <= 0) {
			msg_error("BaseBeamForceField") << "Rectangular dimensions must be positive: "
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

		msg_info("BaseBeamForceField") << "Rectangular cross-section: Ly=" << Ly << ", Lz=" << Lz;
	}

	template<typename DataTypes>
	bool BaseBeamForceField<DataTypes>::validateParameters() const {
		// Check material properties
		if (d_youngModulus.getValue() <= 0) {
			msg_error("BaseBeamForceField") << "Young's modulus must be positive: " << d_youngModulus.getValue();
			return false;
		}

		const Real nu = d_poissonRatio.getValue();
		if (nu <= -1.0 || nu >= 0.5) {
			msg_error("BaseBeamForceField") << "Poisson's ratio must be in range (-1, 0.5): " << nu;
			return false;
		}

		// Check geometric parameters
		const std::string &shape = d_crossSectionShape.getValue().getSelectedItem();

		if (shape == "circular") {
			if (d_radius.getValue() <= 0) {
				msg_error("BaseBeamForceField") << "Radius must be positive: " << d_radius.getValue();
				return false;
			}
		} else if (shape == "hollow_circular") {
			if (d_radius.getValue() <= 0) {
				msg_error("BaseBeamForceField") << "External radius must be positive: " << d_radius.getValue();
				return false;
			}
			if (d_innerRadius.getValue() <= 0) {
				msg_error("BaseBeamForceField") << "Inner radius must be positive for hollow circular section: " << d_innerRadius.getValue();
				return false;
			}
			if (d_innerRadius.getValue() >= d_radius.getValue()) {
				msg_error("BaseBeamForceField") << "Inner radius must be smaller than external radius: "
												<< "rInner=" << d_innerRadius.getValue() << ", r=" << d_radius.getValue();
				return false;
			}
		} else if (shape == "rectangular") {
			if (d_lengthY.getValue() <= 0 || d_lengthZ.getValue() <= 0) {
				msg_error("BaseBeamForceField") << "Rectangular dimensions must be positive: "
												<< "Ly=" << d_lengthY.getValue() << ", Lz=" << d_lengthZ.getValue();
				return false;
			}
		}

		return true;
	}

	template<typename DataTypes>
	bool BaseBeamForceField<DataTypes>::isValidConfiguration() const {
		return validateParameters() && (m_crossSectionArea > 0) && (J > 0) && (Iy > 0) && (Iz > 0);
	}

	template<typename DataTypes>
	void BaseBeamForceField<DataTypes>::printDebugInfo() const {
		msg_info("BaseBeamForceField") << "Cross-section properties:";
		msg_info("BaseBeamForceField") << "  Shape: " << d_crossSectionShape.getValue().getSelectedItem();
		msg_info("BaseBeamForceField") << "  Area: " << m_crossSectionArea;
		msg_info("BaseBeamForceField") << "  J (torsion): " << J;
		msg_info("BaseBeamForceField") << "  Iy (bending): " << Iy;
		msg_info("BaseBeamForceField") << "  Iz (bending): " << Iz;
		msg_info("BaseBeamForceField") << "  E (Young): " << d_youngModulus.getValue();
		msg_info("BaseBeamForceField") << "  ν (Poisson): " << d_poissonRatio.getValue();
		msg_info("BaseBeamForceField") << "  G (Shear): " << getShearModulus();
	}

} // namespace sofa::component::forcefield

