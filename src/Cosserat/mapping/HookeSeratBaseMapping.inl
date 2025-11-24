/*
 * HookeSeratBaseMapping.inl
 * Implementation details for the HookeSeratBaseMapping class.
 * This file contains implementations for functions inline in the HookeSeratBaseMapping class.
 */
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <iostream>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>

namespace Cosserat::mapping {

	using namespace sofa::component::cosserat::liegroups;
	using sofa::helper::getReadAccessor;
	using sofa::type::Quat;
	using sofa::type::Vec3;
	using sofa::type::Vec6;
	using sofa::type::vector;

	template<class TIn1, class TIn2, class TOut>
	HookeSeratBaseMapping<TIn1, TIn2, TOut>::HookeSeratBaseMapping() :
		d_curv_abs_section(initData(&d_curv_abs_section, "curv_abs_input",
									"Curvilinear abscissa of the input sections along the rod")),
		d_curv_abs_frames(initData(&d_curv_abs_frames, "curv_abs_output",
								   "Curvilinear abscissa of the output frames along the rod")),
		d_debug(initData(&d_debug, false, "debug", "Enable debug output")), m_strain_state(nullptr),
		m_rigid_base(nullptr), m_frames(nullptr) {
		msg_info("HookeSeratBaseMapping") << "HookeSeratBaseMapping constructor called  !!!";
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::init() {
		// Initialize pointers to nullptr
		msg_info("HookeSeratBaseMapping") << "Initializing HookeSeratBaseMapping...";

		m_strain_state = nullptr;
		m_rigid_base = nullptr;
		m_frames = nullptr;

		// Check if all required models are present
		if (this->fromModels1.empty()) {
			msg_error() << "Input1 (strain state) not found";
			return;
		}

		if (this->fromModels2.empty()) {
			msg_error() << "Input2 (rigid base) not found";
			return;
		}

		if (this->toModels.empty()) {
			msg_error() << "Output (frames) missing";
			return;
		}

		// Assign mechanical states
		m_strain_state = this->fromModels1[0];
		m_rigid_base = this->fromModels2[0];
		m_frames = this->toModels[0];

		// Get the initial configuration g(X):frames and initialize FrameInfo objects
		if (m_frames) {
			auto xfromData = m_frames->read(sofa::core::vec_id::read_access::position);
			const auto &xfrom = xfromData->getValue();

			// Initialize frame properties using the initial frame states
			const auto frame_count = xfrom.size();

			m_frameProperties.clear();
			m_frameProperties.reserve(frame_count);

			for (size_t i = 0; i < frame_count; ++i) {
				// Convert SOFA coordinates to SE3 transformations
				const auto &frame_i = xfrom[i];
				Vector3 translation(frame_i.getCenter()[0], frame_i.getCenter()[1], frame_i.getCenter()[2]);

				// Convert quaternion to rotation matrix
				const auto &quat = frame_i.getOrientation();
				SO3Type rotation;
				// Convert SOFA quaternion to our SO3 representation
				// SOFA quaternions use [x, y, z, w] order, Eigen uses [w, x, y, z]
				Eigen::Quaternion<double> eigenQuat(quat[3], quat[0], quat[1], quat[2]);
				rotation = SO3Type(eigenQuat.toRotationMatrix());

				SE3Type gXi(rotation, translation);

			// Frame info initialized

				// Create FrameInfo with initial transformation
				// Length and kappa will be set later in initializeFrameProperties
				FrameInfo frameInfo;
				frameInfo.setTransformation(gXi);
				m_frameProperties.push_back(frameInfo);
			}
		}

		// Initialize geometry information
		updateGeometryInfo();

		// Initialize section and frame properties based on geometry
		initializeSectionProperties();

		// No need anymore, this is done in updateGeometryInfo
		initializeFrameProperties();

		// Validation
		if (!validateSectionProperties()) {
			msg_error() << "Invalid section properties detected";
			return;
		}

		// Check continuity (with warning only)
		// if (!checkContinuity()) {
		// 	msg_warning() << "Rod sections are not continuous";
		// }

		// Pre-compute adjoint matrices for performance
		for (const auto &section: m_section_properties) {
			section.getAdjoint(); // Force computation and caching
		}

		// Call parent initialization
		Inherit::init();
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::updateGeometryInfo() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto &curv_abs_frames = d_curv_abs_frames.getValue();

		if (curv_abs_frames.empty()) {
			msg_warning("HookeSeratBaseMapping") << "Empty frames data";
			return;
		}

		const auto frame_count = curv_abs_frames.size();

		// Pré-allocation efficace
		reserveContainers(frame_count);

		size_t current_section_index = 1;
		constexpr double TOLERANCE = 1e-3;

		for (size_t i = 0; i < frame_count; ++i) {
			const auto frame_pos = curv_abs_frames[i];

			// Find the index of the section of on which each frame belong on
			auto result = findSectionIndex(frame_pos, curv_abs_section, current_section_index, TOLERANCE);
			// update the beam's frames strcut
			updateFrameData(i, result.index_for_frame, frame_pos, curv_abs_section);

			// Mettre à jour pour la prochaine itération
			current_section_index = result.index_for_next;
		}
		logCompletionInfo();
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::initializeSectionProperties() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto node_count = curv_abs_section.size();
		const auto &strain = m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

		// Initialize section properties based on geometry
		m_section_properties.clear();
		m_section_properties.reserve(node_count);

		// Save node 0 info
		TangentVector init_strain = TangentVector::Zero(); // Initial curvature, will be updated
		SectionInfo node_0(0., init_strain, 0, SE3Type::computeIdentity()); // First node has zero length
		// The first node is often fix or attached to an object
		m_section_properties.push_back(node_0);

		// Initial strain vector, can be modified later
		TangentVector strain_0 = TangentVector::Zero();

		// compute the length of each beam segment.
		//Initialize a 6D strain vector (angular and linear components)
		std::adjacent_difference(curv_abs_section.begin() + 1, curv_abs_section.end(),
								 std::back_inserter(m_beam_length_vectors));

		for (size_t i = 0; i < node_count - 1; ++i) {
			double length = m_beam_length_vectors[i];
			// Fill with

			SectionInfo node(length, strain_0, i, SE3Type::computeIdentity());
			node.setStrain(strain[i]);
			m_section_properties.push_back(node);
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::initializeFrameProperties() {
		const auto &curv_abs_section = d_curv_abs_section.getValue();
		const auto &curv_abs_frames = d_curv_abs_frames.getValue();

		// Initialize frame properties based on the number of frames
		// The initiation is done in previous updateGeometryInfo

		// Create frame properties for each frame
		for (size_t i = 0; i < curv_abs_frames.size()-1; ++i) {
			if (i < m_indices_vectors.size()) {
				// Calculate frame section length based on position relative to beam nodes
				// This should be the distance from the frame to the closest beam node toward the base
				double frame_length = curv_abs_frames[i + 1] - curv_abs_frames[i];

				// Ensure frameLength is positive (safety check)
				if (frame_length <= 0) {
					msg_warning() << "Frame " << i << " has non-positive length " << frame_length
								  << ". Frame pos: " << curv_abs_frames[i]
								  << ", Section pos: " << curv_abs_section[m_indices_vectors[i] - 1]
								  << ". Using curv_abs_frames[i] instead.";
					frame_length = curv_abs_frames[i];
				}

				m_frameProperties[i+1].setLength(frame_length);
			}
		}
	}

	template<class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::computeTangExpImplementation(const double& curv_abs,
	const TangentVector & strain, const AdjointMatrix &adjoint_matrix, AdjointMatrix & tang_adjoint_matrix)
	{
		// Use Lie group right Jacobian instead of manual trigonometric series
		// This replaces the complex manual computation with the robust Lie group library implementation
		TangentVector scaled_strain = strain * curv_abs;
		tang_adjoint_matrix = SE3Type::rightJacobian(scaled_strain);
	}

	/**
	 * @brief Legacy implementation using manual trigonometric series (kept for verification)
	 * This method computes the tangent exponential map using the original trigonometric series expansion.
	 * Kept for comparison and validation against the new Lie group implementation.
	 */
	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::computeTangExpImplementationLegacy(const double& curv_abs,
		const TangentVector & strain, const AdjointMatrix &adjoint_matrix, AdjointMatrix & tang_adjoint_matrix)
	{
		SReal theta = Vector3(strain(0), strain(1), strain(2)).norm();
		//old method
		//Matrix3 tilde_k = SO3Type::buildAntisymmetric(Vector3(strain_i(0), strain_i(1), strain_i(2)));// getTildeMatrix(Vec3(strain_i(0), strain_i(1), strain_i(2)));
		//Matrix3 tilde_q = SO3Type::buildAntisymmetric(Vector3(strain_i(3), strain_i(4), strain_i(5)));
		//buildAdjoint(tilde_k, tilde_q, ad_Xi);

		//@info: new method with Lie algebra
		//AdjointMatrix adjoint_matrix = node_info.getAdjoint();


		//@todo : compare the result of these two methods for computing
		//@todo : adjoint matrix;

		tang_adjoint_matrix = AdjointMatrix::Zero();
		AdjointMatrix Id6 = AdjointMatrix::Identity();

		if (theta <= std::numeric_limits<double>::epsilon()) {
			double scalar0 = std::pow(curv_abs, 2) / 2.0;
			tang_adjoint_matrix = curv_abs * Id6 + scalar0 * adjoint_matrix;
		} else {
			double scalar1 = (4.0 - 4.0 * cos(curv_abs * theta) -
							  curv_abs * theta * sin(curv_abs * theta)) /
							 (2.0 * theta * theta);
			double scalar2 = (4.0 * curv_abs * theta +
							  curv_abs * theta * cos(curv_abs * theta) -
							  5.0 * sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta);
			double scalar3 = (2.0 - 2.0 * cos(curv_abs * theta) -
							  curv_abs * theta * sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta * theta);
			double scalar4 = (2.0 * curv_abs * theta +
							  curv_abs * theta * cos(curv_abs * theta) -
							  3.0 * sin(curv_abs * theta)) /
							 (2.0 * theta * theta * theta * theta * theta);

			tang_adjoint_matrix = curv_abs * Id6 + scalar1 * adjoint_matrix + scalar2 * adjoint_matrix * adjoint_matrix +
				  scalar3 * adjoint_matrix * adjoint_matrix * adjoint_matrix +
				  scalar4 * adjoint_matrix * adjoint_matrix * adjoint_matrix * adjoint_matrix;
		}
	}

	/**
	 * @brief Test function to compare new and legacy implementations
	 * @return true if implementations produce equivalent results within tolerance
	 */
	template<class TIn1, class TIn2, class TOut>
	bool HookeSeratBaseMapping<TIn1, TIn2, TOut>::testTangExpImplementationEquivalence(const double& curv_abs,
		const TangentVector & strain, const AdjointMatrix &adjoint_matrix, double tolerance = 1e-6)
	{
		AdjointMatrix new_result, legacy_result;

		// Compute using new Lie group implementation
		computeTangExpImplementation(curv_abs, strain, adjoint_matrix, new_result);

		// Compute using legacy trigonometric implementation
		computeTangExpImplementationLegacy(curv_abs, strain, adjoint_matrix, legacy_result);

		// Compare results
		AdjointMatrix diff = new_result - legacy_result;
		double max_diff = diff.cwiseAbs().maxCoeff();

		if (max_diff > tolerance) {
			msg_warning("HookeSeratBaseMapping") << "Tangent exponential implementations differ by " << max_diff
				<< " (tolerance: " << tolerance << ")";
			msg_warning("HookeSeratBaseMapping") << "New result:\n" << new_result;
			msg_warning("HookeSeratBaseMapping") << "Legacy result:\n" << legacy_result;
			return false;
		}

		return true;
	}


	template<class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::updateTangExpSE3() {

		auto node_count = m_section_properties.size();

		//update node's tang SE3 matrix
		AdjointMatrix tang_matrix = AdjointMatrix::Zero();
		m_section_properties[0].setTanAdjointMatrix(tang_matrix);

		for (auto i = 1; i< node_count; i++ ) {
			auto node_info = m_section_properties[i];
			computeTangExpImplementation(
				node_info.getLength(),
				node_info.getStrainsVec(),
				node_info.getAdjoint(),
				tang_matrix);
			node_info.setTanAdjointMatrix(tang_matrix);
			std::cout << "Node[" << i << "] tang adjoint matrix: \n" << tang_matrix << std::endl;
		}

		//update frames's tang SE3 matrix
		auto frame_count = m_frameProperties.size();
		for (auto i = 0; i<frame_count; i++) {
			auto frame_info = m_frameProperties[i];
			auto related_section_index = m_frameProperties[i].get_related_beam_index_();
			auto frame_strain = m_section_properties[related_section_index].getStrainsVec();
			computeTangExpImplementation(
				frame_info.getDistanceToNearestBeamNode(),
				frame_strain,
				frame_info.getAdjoint(),
				tang_matrix);
			frame_info.setTanAdjointMatrix(tang_matrix);
			std::cout << "Frame[" << i << "] tang adjoint matrix: \n"
			<< tang_matrix << std::endl;
		}
	}

	// Implementation of BeamTopology
	bool BeamTopology::isValid() const {
		if (parent_indices.empty()) return true;

		// Check for valid tree structure (no cycles, valid parent indices)
		std::vector<bool> visited(parent_indices.size(), false);
		std::vector<int> parent_count(parent_indices.size(), 0);

		// Count parents
		for (int parent : parent_indices) {
			if (parent >= 0 && static_cast<size_t>(parent) < parent_indices.size()) {
				parent_count[parent]++;
			} else if (parent != -1) {
				return false; // Invalid parent index
			}
		}

		// Check for cycles using DFS
		std::function<bool(size_t)> hasCycle = [&](size_t node) -> bool {
			if (visited[node]) return true; // Cycle detected
			visited[node] = true;

			for (size_t child : getChildren(node)) {
				if (hasCycle(child)) return true;
			}

			visited[node] = false; // Backtrack
			return false;
		};

		// Check each node for cycles
		std::fill(visited.begin(), visited.end(), false);
		for (size_t i = 0; i < parent_indices.size(); ++i) {
			if (hasCycle(i)) return false;
		}

		return true;
	}

	std::vector<size_t> BeamTopology::getChildren(size_t section_idx) const {
		std::vector<size_t> children;
		for (size_t i = 0; i < parent_indices.size(); ++i) {
			if (static_cast<size_t>(parent_indices[i]) == section_idx) {
				children.push_back(i);
			}
		}
		return children;
	}

	// Implementation of BeamStateEstimator
	BeamStateEstimator::BeamStateEstimator()
		: process_noise_(Eigen::Matrix<double, 12, 12>::Identity() * 0.01),
		  measurement_noise_(Eigen::Matrix<double, 6, 6>::Identity() * 0.1) {
		// Initialize with identity pose and zero strain
		SE3Type identity_pose = SE3Type::computeIdentity();
		StrainState zero_strain(SO3Type(Vector3::Zero()),
							   sofa::component::cosserat::liegroups::RealSpace<double, 3>(Vector3::Zero()));

		Eigen::Matrix<double, 12, 12> initial_cov = Eigen::Matrix<double, 12, 12>::Identity() * 0.1;
		initialize(identity_pose, zero_strain, initial_cov);
	}

	void BeamStateEstimator::initialize(const SE3Type& initial_pose, const StrainState& initial_strain,
									   const Eigen::Matrix<double, 12, 12>& initial_covariance) {
		pose_estimate_ = sofa::component::cosserat::liegroups::GaussianOnManifold<SE3Type>(initial_pose,
			initial_covariance.template topLeftCorner<6, 6>());
		strain_estimate_ = initial_strain;
		state_covariance_ = initial_covariance;
		initialized_ = true;
	}

	void BeamStateEstimator::predict(const TangentVector& control_input, double dt) {
		if (!initialized_) return;

		// Simple prediction model: constant velocity with process noise
		// In a real implementation, this would use the beam dynamics model
		state_covariance_ += process_noise_ * dt;
		// Pose prediction would use the control input through the beam model
		// For now, we keep the pose estimate unchanged
	}

	void BeamStateEstimator::update(const SE3Type& measurement,
								   const Eigen::Matrix<double, 6, 6>& measurement_covariance) {
		if (!initialized_) return;

		// Kalman update for pose measurement
		auto pose_cov = state_covariance_.template topLeftCorner<6, 6>();
		Eigen::Matrix<double, 6, 6> innovation_cov = pose_cov + measurement_covariance;

		// Compute Kalman gain
		Eigen::Matrix<double, 6, 6> kalman_gain = pose_cov * innovation_cov.inverse();

		// Update pose estimate (simplified - would need proper Lie group update)
		// This is a placeholder for the actual Lie group Kalman update
		auto pose_cov_updated = (Eigen::Matrix<double, 6, 6>::Identity() - kalman_gain) * pose_cov;
		state_covariance_.template topLeftCorner<6, 6>() = pose_cov_updated;
	}

	void BeamStateEstimator::updateStrain(const StrainState& strain_measurement,
										 const Eigen::Matrix<double, 6, 6>& measurement_covariance) {
		if (!initialized_) return;

		// Kalman update for strain measurement
		auto strain_cov = state_covariance_.template bottomRightCorner<6, 6>();
		Eigen::Matrix<double, 6, 6> innovation_cov = strain_cov + measurement_covariance;

		// Compute Kalman gain
		Eigen::Matrix<double, 6, 6> kalman_gain = strain_cov * innovation_cov.inverse();

		// Update strain estimate (simplified)
		auto strain_cov_updated = (Eigen::Matrix<double, 6, 6>::Identity() - kalman_gain) * strain_cov;
		state_covariance_.template bottomRightCorner<6, 6>() = strain_cov_updated;
	}

	double BeamStateEstimator::getEstimationConfidence() const {
		if (!initialized_) return 0.0;
		return 1.0 / state_covariance_.trace();
	}

	void BeamStateEstimator::reset() {
		initialized_ = false;
		state_covariance_.setZero();
	}

	// HookeSeratBaseMapping multi-section and state estimation methods
	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::setBeamTopology(const BeamTopology& topology) {
		if (!topology.isValid()) {
			msg_error() << "Invalid beam topology provided";
			return;
		}
		m_beam_topology = topology;
		m_multi_section_enabled = true;
		msg_info() << "Beam topology set with " << topology.getNumSections() << " sections";
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::enableStateEstimation(bool enable) {
		if (enable && !m_state_estimator) {
			m_state_estimator = std::make_unique<BeamStateEstimator>();
			msg_info() << "State estimation enabled";
		} else if (!enable && m_state_estimator) {
			m_state_estimator.reset();
			msg_info() << "State estimation disabled";
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::updateStateEstimate(const SE3Type& measurement,
																	 const Eigen::Matrix<double, 6, 6>& covariance) {
		if (m_state_estimator) {
			m_state_estimator->update(measurement, covariance);
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::predictState(const TangentVector& control_input, double dt) {
		if (m_state_estimator) {
			m_state_estimator->predict(control_input, dt);
		}
	}

	template<class TIn1, class TIn2, class TOut>
	double HookeSeratBaseMapping<TIn1, TIn2, TOut>::getEstimationConfidence() const {
		return m_state_estimator ? m_state_estimator->getEstimationConfidence() : 0.0;
	}

	// Performance optimization methods
	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::enableParallelComputation(bool enable) {
		parallel_computation_enabled_ = enable;
		if (enable) {
			// Determine optimal thread count based on system capabilities
			optimal_thread_count_ = std::max(1u, std::thread::hardware_concurrency() / 2);
			msg_info() << "Parallel computation enabled with " << optimal_thread_count_ << " threads";
		} else {
			optimal_thread_count_ = 1;
			msg_info() << "Parallel computation disabled";
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::clearComputationCache() {
		computation_cache_.clear();
		last_cache_clear_ = std::chrono::steady_clock::now();
		msg_info() << "Computation cache cleared";
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::runPerformanceBenchmark(size_t iterations) {
		msg_info() << "Running performance benchmark with " << iterations << " iterations";
		m_jacobian_stats.benchmarkJacobianComputation(iterations);
		msg_info() << "Performance benchmark completed";
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::printPerformanceReport() const {
		m_jacobian_stats.printPerformanceReport();
	}

	// JacobianStats benchmarking implementation
	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::JacobianStats::benchmarkJacobianComputation(size_t iterations) {
		reset();

		// Create test data
		TangentVector test_strain = TangentVector::Random();

		for (size_t i = 0; i < iterations; ++i) {
			startTiming();
			AdjointMatrix jacobian = SE3Type::rightJacobian(test_strain);
			endTiming();

			// Use the result to prevent optimization
			jacobian(0, 0) += 0.0;
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratBaseMapping<TIn1, TIn2, TOut>::JacobianStats::printPerformanceReport() const {
		std::cout << "=== Jacobian Computation Performance Report ===" << std::endl;
		std::cout << "Total computations: " << computation_count << std::endl;
		std::cout << "Cache hits: " << cache_hits << std::endl;
		std::cout << "Cache hit rate: " << cacheHitRate() * 100.0 << "%" << std::endl;
		std::cout << "Average computation time: " << averageComputationTime() << " ms" << std::endl;
		std::cout << "Total computation time: " << total_computation_time << " ms" << std::endl;
		std::cout << "Cache size: " << jacobian_cache.size() << " entries" << std::endl;
		std::cout << "===============================================" << std::endl;
	}

} // namespace Cosserat::mapping
