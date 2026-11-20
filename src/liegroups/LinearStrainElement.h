#pragma once

/**
 * @file LinearStrainElement.h
 * @brief Linear Strain Element (LSE) of a Cosserat rod: nodal poses on SE(3) + one strain slope.
 *
 * Reference: L. Xun, B. Tamadazte, "Geometrically explicit Cosserat-rod modeling with piecewise
 * linear strain for complex rod systems", Computers & Structures 325 (2026) 108181.
 * Equation numbers below refer to that paper.
 *
 * ## Model
 * An element of length h joins two nodal poses g_a, g_b in SE(3). Its strain field (body frame,
 * g' = g ξ̂, ξ = [κ ; ε] angular first) varies linearly:
 *
 *     ξ(s) = ξ̄ + (s − h/2) β ,   s in [0, h]                                   (3)
 *
 * with ξ̄ the mean strain and β the strain slope. Integrating g' = g ξ̂ with the fourth-order
 * Magnus expansion is exact in closed form for a linear strain:
 *
 *     g_b = g_a Exp(Ω),   Ω = A ξ̄,   A = h·I − (h³/12) ad_β                    (4)(5)
 *
 * so the mean strain is reconstructed *explicitly*, with no iteration:
 *
 *     ξ̄ = A⁻¹ Log(g_a⁻¹ g_b)                                                   (6)
 *
 * (β = 0 gives A = h·I, i.e. the constant-strain element ξ̄ = Log(g_a⁻¹ g_b)/h.)
 *
 * ## Conventions
 *  - ξ is the *full* se(3) strain (a straight, unstretched rod has ε = (1,0,0)). Rest strains
 *    and the `+e_x` offset of the plugin's Vec6 strain are the caller's business; only β (a
 *    difference) is offset-free. ad_ξ̄ in J₃ needs the full ξ̄.
 *  - Perturbations are right/body ones: g ← g·Exp(ζ̂). Hence the Jacobians below are expressed
 *    in the body frame of each node.
 *  - Element degrees of freedom are ordered [ζ_a ; ζ_b ; β] (18 values).
 *
 * ## Jacobians  δξ̄ = J₁ δζ_a + J₂ δζ_b + J₃ δβ
 *     J₂ =  A⁻¹ Jr⁻¹(Ω)
 *     J₁ = −A⁻¹ Jr⁻¹(Ω) Ad_{Exp(−Ω)}
 *     J₃ = −(h³/12) A⁻¹ ad_ξ̄
 * Eq. (8) of the paper prints J₁ without the minus sign; the sign above is the one that matches
 * central finite differences (see test_LinearStrainElement.cpp) and the usual
 * Log(g_a⁻¹ g_b) derivative.
 *
 * ## Energy, residual and Gauss–Newton tangent
 *     U = ½ h (ξ̄−ξ₀)ᵀK(ξ̄−ξ₀) + (h³/24) βᵀKβ                                  (21)
 *     R₁ = h J₁ᵀK(ξ̄−ξ₀),  R₂ = h J₂ᵀK(ξ̄−ξ₀),  R₃ = h J₃ᵀK(ξ̄−ξ₀) + (h³/12) Kβ  (23-25)
 *     K^GN_ij = h J_iᵀ K J_j  for (i,j) ≠ (3,3),  K^GN_33 = h J₃ᵀKJ₃ + (h³/12) K  (30)
 * R is the gradient of U (the internal force is −R). The tangent holds the Jacobians constant
 * (Gauss–Newton): it is exact when ξ̄ = ξ₀ and otherwise drops the term ∝ (ξ̄−ξ₀).
 *
 * The header only needs Eigen (no SOFA), so it can be unit-tested stand-alone.
 */

#include <liegroups/SE3.h>

#include <Eigen/Dense>

namespace sofa::component::cosserat::liegroups {

template<typename _Scalar = double>
class LinearStrainElement {
public:
	using Scalar = _Scalar;
	using SE3Type = SE3<Scalar>;
	using Vec6 = typename SE3Type::TangentVector;
	using Mat6 = typename SE3Type::AdjointMatrix;
	using Vec18 = Eigen::Matrix<Scalar, 18, 1>;
	using Mat18 = Eigen::Matrix<Scalar, 18, 18>;

	/// Everything the mapping needs about one element at the current configuration.
	struct Kinematics {
		Vec6 omega; ///< Log(g_a⁻¹ g_b): twist integrated over the element
		Mat6 A;     ///< h·I − (h³/12) ad_β
		Vec6 xiBar; ///< mean strain, Eq. (6)
		Mat6 J1;    ///< ∂ξ̄/∂ζ_a (body frame of g_a)
		Mat6 J2;    ///< ∂ξ̄/∂ζ_b (body frame of g_b)
		Mat6 J3;    ///< ∂ξ̄/∂β
	};

	/// Gradient of the element energy with respect to (ζ_a, ζ_b, β).
	struct Residual {
		Vec6 R1;
		Vec6 R2;
		Vec6 R3;
		Vec18 stacked() const {
			Vec18 r;
			r << R1, R2, R3;
			return r;
		}
	};

	/// A = h·I − (h³/12) ad_β. Always invertible: ad_β has the eigenvalues of [κ_β]× (0, ±i‖κ_β‖),
	/// hence det A = h⁶ (1 + (h²‖κ_β‖/12)²)² > 0.
	static Mat6 computeA(const Vec6 &beta, Scalar h) {
		return h * Mat6::Identity() - (h * h * h / Scalar(12)) * SE3Type::computeAd(beta);
	}

	/// Ω = A ξ̄, Eq. (5).
	static Vec6 integratedTwist(const Vec6 &xiBar, const Vec6 &beta, Scalar h) {
		return computeA(beta, h) * xiBar;
	}

	/// g_b = g_a Exp(A ξ̄), Eq. (4)-(5): the far node of an element given its strain field.
	static SE3Type endPose(const SE3Type &ga, const Vec6 &xiBar, const Vec6 &beta, Scalar h) {
		return ga.computeCompose(SE3Type::computeExp(integratedTwist(xiBar, beta, h)));
	}

	/// ξ̄ = A⁻¹ Log(g_a⁻¹ g_b), Eq. (6). With β = 0 this is Log(g_a⁻¹ g_b)/h, which is how a rest
	/// strain ξ₀ is computed from the initial nodal poses.
	static Vec6 meanStrain(const SE3Type &ga, const SE3Type &gb, const Vec6 &beta, Scalar h) {
		const Vec6 omega = relativeLog(ga, gb);
		return computeA(beta, h).partialPivLu().solve(omega);
	}

	/// Mean strain, Jacobians J₁, J₂, J₃ (Eq. 7-10) and intermediate quantities.
	static Kinematics evaluate(const SE3Type &ga, const SE3Type &gb, const Vec6 &beta, Scalar h) {
		Kinematics k;
		k.omega = relativeLog(ga, gb);
		k.A = computeA(beta, h);
		const Eigen::PartialPivLU<Mat6> lu(k.A);
		k.xiBar = lu.solve(k.omega);

		const Mat6 JrInv = SE3Type::computeRightJacobianInverse(k.omega);
		const Mat6 AdExpMinusOmega = SE3Type::computeExp(-k.omega).computeAdjoint();
		const Mat6 AinvJrInv = lu.solve(JrInv);

		k.J2 = AinvJrInv;
		k.J1 = -AinvJrInv * AdExpMinusOmega;
		k.J3 = -(h * h * h / Scalar(12)) * lu.solve(Mat6(SE3Type::computeAd(k.xiBar)));
		return k;
	}

	/// Element energy, Eq. (21).
	static Scalar energy(const Vec6 &xiBar, const Vec6 &xi0, const Vec6 &beta, const Mat6 &K, Scalar h) {
		const Vec6 d = xiBar - xi0;
		return Scalar(0.5) * h * d.dot(K * d) + h * h * h / Scalar(24) * beta.dot(K * beta);
	}

	/// Gradient of the element energy, Eq. (23)-(25).
	static Residual residual(const Kinematics &k, const Vec6 &xi0, const Vec6 &beta, const Mat6 &K, Scalar h) {
		const Vec6 Kd = K * (k.xiBar - xi0);
		Residual r;
		r.R1 = h * k.J1.transpose() * Kd;
		r.R2 = h * k.J2.transpose() * Kd;
		r.R3 = h * k.J3.transpose() * Kd + h * h * h / Scalar(12) * (K * beta);
		return r;
	}

	/// Gauss–Newton tangent stiffness, Eq. (30), ordered [ζ_a ; ζ_b ; β].
	static Mat18 gaussNewtonTangent(const Kinematics &k, const Mat6 &K, Scalar h) {
		const Mat6 *J[3] = {&k.J1, &k.J2, &k.J3};
		Mat18 T = Mat18::Zero();
		for (int i = 0; i < 3; ++i) {
			for (int j = 0; j < 3; ++j) {
				T.template block<6, 6>(6 * i, 6 * j) = h * J[i]->transpose() * K * (*J[j]);
			}
		}
		T.template block<6, 6>(12, 12) += h * h * h / Scalar(12) * K;
		return T;
	}

private:
	static Vec6 relativeLog(const SE3Type &ga, const SE3Type &gb) {
		return ga.computeInverse().computeCompose(gb).computeLog();
	}
};

using LinearStrainElementd = LinearStrainElement<double>;

} // namespace sofa::component::cosserat::liegroups
