#pragma once

#include <Cosserat/config.h>
#include <liegroups/SE3.h>
#include <vector>

namespace Cosserat::mapping {

	using namespace sofa::component::cosserat::liegroups;

	/**
	 * @brief Class for interpolating between two beam shapes (configurations).
	 *
	 * This class provides methods to interpolate between two sets of SE(3) frames
	 * representing the shape of a beam, using geodesic interpolation on the SE(3) manifold.
	 */
	class BeamShapeInterpolation {
	public:
		using SE3Type = SE3<double>;

		/**
		 * @brief Interpolate between two beam shapes.
		 *
		 * @param shape1 First shape (vector of SE3 frames).
		 * @param shape2 Second shape (vector of SE3 frames).
		 * @param t Interpolation parameter (0.0 to 1.0).
		 * @return Interpolated shape.
		 */
		static std::vector<SE3Type> interpolateShapes(const std::vector<SE3Type> &shape1,
													  const std::vector<SE3Type> &shape2, double t) {
			if (shape1.size() != shape2.size()) {
				// Handle error or return empty
				return {};
			}

			std::vector<SE3Type> result;
			result.reserve(shape1.size());

			for (size_t i = 0; i < shape1.size(); ++i) {
				// Geodesic interpolation: A * Exp(t * Log(A^-1 * B))
				result.push_back(shape1[i].interpolate(shape2[i], t));
			}

			return result;
		}

		/**
		 * @brief Blend multiple shapes with weights.
		 *
		 * @param shapes List of shapes.
		 * @param weights List of weights (should sum to 1).
		 * @return Blended shape.
		 */
		static std::vector<SE3Type> blendShapes(const std::vector<std::vector<SE3Type>> &shapes,
												const std::vector<double> &weights) {
			if (shapes.empty() || shapes.size() != weights.size()) {
				return {};
			}

			size_t num_frames = shapes[0].size();
			for (const auto &shape: shapes) {
				if (shape.size() != num_frames)
					return {};
			}

			std::vector<SE3Type> result;
			result.reserve(num_frames);

			// For each frame index
			for (size_t i = 0; i < num_frames; ++i) {
				// Compute weighted average on manifold (Fréchet mean)
				// For SE(3), a simple approximation is iterative or just linear blending in tangent space of the first
				// shape. Let's use tangent space blending relative to the first shape (shapes[0]).

				SE3Type reference = shapes[0][i];
				SE3Type::TangentVector weighted_sum = SE3Type::TangentVector::Zero();
				double total_weight = 0.0;

				for (size_t k = 0; k < shapes.size(); ++k) {
					// Log map relative to reference
					// v_k = Log(ref^-1 * shape_k)
					SE3Type diff = reference.computeInverse() * shapes[k][i];
					weighted_sum += weights[k] * diff.log();
					total_weight += weights[k];
				}

				if (std::abs(total_weight) > 1e-6) {
					weighted_sum /= total_weight;
				}

				// Result = ref * Exp(weighted_sum)
				result.push_back(reference * SE3Type::computeExp(weighted_sum));
			}

			return result;
		}
	};

} // namespace Cosserat::mapping
