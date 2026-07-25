/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                                                                             *
 * This file is part of the Cosserat plugin.                                   *
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <chrono>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>
#include <liegroups/Types.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

#include <numeric>
#include <stdexcept>
#include <vector>

namespace Cosserat::mapping {

	// ──────────────────────────────────────────────────────────────────────────
	// Type aliases shared across all mapping classes (moved here from the
	// formerly monolithic CosseratGeometryMapping.h to keep them close to the
	// geometry helpers that consume them).
	// ──────────────────────────────────────────────────────────────────────────
	using SE3Type        = sofa::component::cosserat::liegroups::SE3<double>;
	using SO3Type        = sofa::component::cosserat::liegroups::SO3<double>;
	using Vector3        = typename SE3Type::Vector3;
	using Matrix3        = typename SE3Type::Matrix3;
	using Matrix4        = typename SE3Type::Matrix4;
	using AdjointMatrix  = typename SE3Type::AdjointMatrix;
	using JacobianMatrix = typename SE3Type::JacobianMatrix;
	using TangentVector  = typename SE3Type::TangentVector;

	// ──────────────────────────────────────────────────────────────────────────
	// SectionInfo / FrameInfo / BeamTopology
	//
	// These describe a single discretised section, a single output frame, and
	// the parent–child topology of a multi-section beam respectively. They are
	// SOFA-agnostic value types (no BaseObject inheritance).
	// ──────────────────────────────────────────────────────────────────────────

	/**
	 * @brief Class encapsulating the properties of a Cosserat beam node
	 *
	 * STRAIN CONVENTION:
	 * Strains follow Cosserat convention: [φx, φy, φz, ρx, ρy, ρz]ᵀ
	 *   - angular_strain_ = φ = [φx, φy, φz] (head<3>): torsion, bending Y, bending Z
	 *   - linear_strain_  = ρ = [ρx, ρy, ρz] (tail<3>): elongation, shearing Y, shearing Z
	 * See liegroups/STRAIN_CONVENTION.md for detailed documentation.
	 */
	class SectionInfo {
	private:
		double sec_length_ = 0.0;
		Vector3 angular_strain_ = Vector3::Zero(); ///< φ = [φx, φy, φz]: torsion, bending Y, bending Z
		Vector3 linear_strain_  = Vector3::Zero(); ///< ρ = [ρx, ρy, ρz]: elongation, shearing Y, shearing Z
		TangentVector strain_   = TangentVector::Zero(); ///< Full strain [φ, ρ]ᵀ ∈ se(3)

		unsigned int index_0_ = 0;
		unsigned int index_1_ = index_0_ + 1;

		SE3Type gX_;

		// Lazy cache of computeAdjoint() / computeCoAdjoint(). Invalidated by setTransformation().
		mutable AdjointMatrix adjoint_;
		mutable AdjointMatrix coAdjoint_;
		mutable bool adjoint_dirty_   = true;
		mutable bool coAdjoint_dirty_ = true;
		AdjointMatrix tang_adjoint_;
		AdjointMatrix tang_co_adjoint_;

	public:
		SectionInfo() = default;

		SectionInfo(double length, const TangentVector &strain, const unsigned int i0,
					const SE3Type &transform = SE3Type::computeIdentity()) :
			sec_length_(length), strain_(strain), index_0_(i0), gX_(transform) {}

		SectionInfo(double length, const Vector3 &angular_strain, const unsigned int i0) :
			sec_length_(length), angular_strain_(angular_strain), index_0_(i0) {}

		double getLength() const { return sec_length_; }
		void setLength(double length) {
			if (length <= 0)
				throw std::invalid_argument("Section length must be positive");
			sec_length_ = length;
		}

		/**
		 * @brief Set strain values from various vector types
		 *
		 * STRAIN INDEXING:
		 * For Vec<6>: strain = [φx, φy, φz, ρx, ρy, ρz]
		 *   - strain[0-2] → angular_strain_ (φ: torsion, bending Y, bending Z)
		 *   - strain[3-5] → linear_strain_ (ρ: elongation, shearing Y, shearing Z)
		 */
		template<typename VecType>
		void setStrain(const VecType &strain) {
			if constexpr (std::is_same_v<VecType, sofa::type::Vec<3, double>>) {
				for (int i = 0; i < 3; ++i) {
					angular_strain_[i] = strain[i];
				}
				strain_.head<3>() = angular_strain_;
				strain_.tail<3>() = linear_strain_;
			} else if constexpr (std::is_same_v<VecType, sofa::type::Vec<6, double>>) {
				for (int i = 0; i < 3; ++i) {
					angular_strain_[i] = strain[i];
					linear_strain_[i]  = strain[i + 3];
				}
				strain_.head<3>() = angular_strain_;
				strain_.tail<3>() = linear_strain_;
			} else {
				if constexpr (requires { VecType::SizeAtCompileTime; }) {
					if constexpr (VecType::SizeAtCompileTime == 3) {
						angular_strain_ = strain;
					} else if constexpr (VecType::SizeAtCompileTime == 6) {
						angular_strain_ = strain.template head<3>();
						linear_strain_  = strain.template tail<3>();
					}
					strain_.head<3>() = angular_strain_;
					strain_.tail<3>() = linear_strain_;
				} else {
					if (strain.size() == 3) {
						for (int i = 0; i < 3; ++i) {
							angular_strain_[i] = strain[i];
						}
					} else if (strain.size() == 6) {
						for (int i = 0; i < 3; ++i) {
							angular_strain_[i] = strain[i];
							linear_strain_[i]  = strain[i + 3];
						}
					}
					strain_.head<3>() = angular_strain_;
					strain_.tail<3>() = linear_strain_;
				}
			}
		}

		const TangentVector &getStrainsVec() const { return strain_; }

		unsigned int getIndex0() const { return index_0_; }
		unsigned int getIndex1() const { return index_1_; }
		void setIndices(unsigned int i0) { index_0_ = i0; index_1_ = i0 + 1; }

		const SE3Type &getTransformation() const { return gX_; }
		void setTransformation(const SE3Type &transform) {
			gX_ = transform;
			adjoint_dirty_   = true;
			coAdjoint_dirty_ = true;
		}

		Matrix4 getTransformationMatrix() const { return gX_.matrix(); }

		void setTransformationFromMatrix(const Matrix4 &matrix) {
			gX_ = SE3Type(matrix);
			adjoint_dirty_   = true;
			coAdjoint_dirty_ = true;
		}

		const AdjointMatrix &getAdjoint() const {
			if (adjoint_dirty_) {
				adjoint_       = gX_.computeAdjoint();
				adjoint_dirty_ = false;
			}
			return adjoint_;
		}

		const AdjointMatrix &getCoAdjoint() const {
			if (coAdjoint_dirty_) {
				coAdjoint_       = gX_.computeCoAdjoint();
				coAdjoint_dirty_ = false;
			}
			return coAdjoint_;
		}

		const AdjointMatrix &getTangAdjointMatrix() const { return tang_adjoint_; }
		void setTanAdjointMatrix(const AdjointMatrix &tang_adjoint_mat) { tang_adjoint_ = tang_adjoint_mat; }

		SE3Type getLocalTransformation(double t) const {
			if (t < 0.0 || t > 1.0)
				throw std::invalid_argument("Parameter t must be in [0,1]");

			TangentVector xi;
			xi.template head<3>() = t * sec_length_ * Vector3::UnitX();
			xi.template tail<3>() = t * angular_strain_;

			return gX_ * SE3Type::computeExp(xi);
		}

		TangentVector getTransformationDerivative(double /*t*/) const {
			TangentVector xi;
			xi.template head<3>() = sec_length_ * Vector3::UnitX();
			xi.template tail<3>() = angular_strain_;
			return getAdjoint() * xi;
		}

		double distanceTo(const SectionInfo &other, double w_rot = 1.0, double w_trans = 1.0) const {
			return gX_.distance(other.gX_, w_rot, w_trans);
		}

		Vector3 transformPoint(const Vector3 &point) const { return gX_.computeAction(point); }

		bool isApprox(const SectionInfo &other, double eps = 1e-6) const {
			return gX_.computeIsApprox(other.gX_, eps) &&
				   (angular_strain_ - other.angular_strain_).norm() < eps &&
				   std::abs(sec_length_ - other.sec_length_) < eps;
		}

		SectionInfo inverse() const {
			return SectionInfo(sec_length_, -strain_, index_0_, gX_.computeInverse());
		}

		SectionInfo compose(const SectionInfo &other) const {
			SE3Type composed_transform = gX_.compose(other.gX_);
			TangentVector composed_strain;
			composed_strain.head<3>() = angular_strain_ + other.angular_strain_;
			composed_strain.tail<3>() = linear_strain_ + other.linear_strain_;
			double total_length = sec_length_ + other.sec_length_;
			return SectionInfo(total_length, composed_strain, index_0_, composed_transform);
		}
	};

	/**
	 * @brief Class for output frame properties along a Cosserat beam.
	 *
	 * Each frame is attached to a parent section (by index) and stored with
	 * its arc-length offset from the section start. Renames applied here vs.
	 * the legacy CosseratGeometryMapping.h:
	 *   related_beam_index_           → related_section_index_
	 *   get_related_beam_index_()     → getRelatedSectionIndex()
	 *   set_related_beam_index_()     → setRelatedSectionIndex()
	 *   distance_to_nearest_beam_node → distance_to_section_start
	 *   getDistanceToNearestBeamNode() → getDistanceToSectionStart()
	 *   setDistanceToNearestBeamNode() → setDistanceToSectionStart()
	 */
	class FrameInfo {
	private:
		double frames_sect_length_ = 0.0;
		TangentVector kappa_       = TangentVector::Zero(); ///< angular strain + linear strain (full 6D)
		unsigned int index_0_      = 0;
		unsigned int index_1_      = 1;
		unsigned int related_section_index_ = 0; ///< Index of the parent section
		double distance_to_section_start    = 0.0; ///< Arc-length offset from the parent section start
		SE3Type transformation_;

		mutable AdjointMatrix adjoint_;
		mutable AdjointMatrix coAdjoint_;
		mutable bool adjoint_dirty_   = true;
		mutable bool coAdjoint_dirty_ = true;
		AdjointMatrix tang_adjoint_;

	public:
		FrameInfo() = default;

		double getLength() const { return frames_sect_length_; }
		void setLength(double length) {
			if (length <= 0)
				throw std::invalid_argument("Frame length must be positive");
			frames_sect_length_ = length;
		}

		unsigned int getRelatedSectionIndex() const { return related_section_index_; }
		void setRelatedSectionIndex(unsigned int index) { related_section_index_ = index; }

		double getDistanceToSectionStart() const { return distance_to_section_start; }
		void setDistanceToSectionStart(double distance) {
			if (distance < 0)
				throw std::invalid_argument("Distance to section start must be non-negative");
			distance_to_section_start = distance;
		}

		const TangentVector &getKappa() const { return kappa_; }
		void setKappa(const TangentVector &k) { kappa_ = k; }

		const SE3Type &getTransformation() const { return transformation_; }
		SE3Type getInverseTransformation() const { return transformation_.inverse(); }

		void setTransformation(const SE3Type &transform) {
			transformation_  = transform;
			adjoint_dirty_   = true;
			coAdjoint_dirty_ = true;
		}

		const AdjointMatrix &getAdjoint() const {
			if (adjoint_dirty_) {
				adjoint_       = transformation_.computeAdjoint();
				adjoint_dirty_ = false;
			}
			return adjoint_;
		}

		const AdjointMatrix &getCoAdjoint() const {
			if (coAdjoint_dirty_) {
				coAdjoint_       = transformation_.computeCoAdjoint();
				coAdjoint_dirty_ = false;
			}
			return coAdjoint_;
		}

		const AdjointMatrix &getTangAdjointMatrix() const { return tang_adjoint_; }
		void setTanAdjointMatrix(const AdjointMatrix &tang_adjoint_mat) { tang_adjoint_ = tang_adjoint_mat; }

		SE3Type getLocalTransformation(double t) const {
			if (t < 0.0 || t > 1.0)
				throw std::invalid_argument("Parameter t must be in [0,1]");

			TangentVector xi = t * kappa_;
			xi.template head<3>() += t * frames_sect_length_ * Vector3::UnitX();
			return transformation_ * SE3Type::computeExp(xi);
		}

		friend std::ostream &operator<<(std::ostream &os, const FrameInfo &frame) {
			os << "FrameInfo{length=" << frame.frames_sect_length_
			   << ", related_section=" << frame.related_section_index_
			   << ", distance_to_section_start=" << frame.distance_to_section_start
			   << ", kappa=[" << frame.kappa_.transpose() << "]"
			   << ", transformation=" << frame.transformation_ << "}";
			return os;
		}
	};

	/**
	 * @brief Tree topology of a multi-section beam.
	 *
	 * Each section may have a parent (via index) and a relative transform from
	 * the parent's frame. Used for branching beams.
	 */
	struct BeamTopology {
		std::vector<int>     parent_indices;
		std::vector<SE3Type> relative_transforms;
		std::vector<double>  connection_stiffnesses;

		bool isValid() const {
			return !parent_indices.empty() && parent_indices.size() == relative_transforms.size();
		}

		std::vector<size_t> getChildren(size_t section_idx) const {
			std::vector<size_t> children;
			for (size_t i = 0; i < parent_indices.size(); ++i) {
				if (parent_indices[i] == static_cast<int>(section_idx)) {
					children.push_back(i);
				}
			}
			return children;
		}

		size_t getNumSections() const { return parent_indices.size(); }
	};

	// ──────────────────────────────────────────────────────────────────────────
	// CosseratBeamGeometry<TStrain>
	//
	// Non-BaseObject mixin holding the geometry-only state and helpers shared
	// by every Cosserat mapping (Strain2Frames, Frames2Strain, …).
	//
	// Inheritance model (Z1):
	//   - This class does NOT inherit from BaseObject (avoids the diamond
	//     when a derived mapping also inherits from Multi2Mapping or Mapping).
	//   - Data fields (d_curv_abs_section / d_curv_abs_frames / d_debug) live
	//     here as members but are REGISTERED via initData() by the derived
	//     class's constructor — the derived class IS a BaseObject and provides
	//     the SOFA-side wiring.
	//   - All logging uses the msg_*("CosseratBeamGeometry") form to avoid
	//     needing `this` to be a BaseObject.
	//
	// Renames applied in this extraction (vs. legacy CosseratGeometryMapping.h):
	//   m_frameProperties        → m_frame_properties
	//   m_indices_vectors        → m_frame_to_section_indices
	//   m_indices_vectors_draw   → m_frame_to_section_indices_draw
	//   m_beam_length_vectors    → m_section_length_vectors
	//   m_frames_length_vectors  → m_frame_distance_to_section_start
	// ──────────────────────────────────────────────────────────────────────────
	template<class TStrain>
	class CosseratBeamGeometry {
	public:
		using Strain = TStrain;
		using StrainCoord = sofa::Coord_t<TStrain>;
		using StrainVecCoord = sofa::VecCoord_t<TStrain>;

		static constexpr bool ENABLE_GEOMETRY_LOGGING = true;

	protected:
		// ── Geometry state ────────────────────────────────────────────────
		std::vector<SectionInfo> m_section_properties;
		std::vector<FrameInfo>   m_frame_properties;
		BeamTopology             m_topology;

		// Frame ↔ section index maps and per-frame arc-length offsets.
		std::vector<unsigned int> m_frame_to_section_indices;
		std::vector<unsigned int> m_frame_to_section_indices_draw;
		std::vector<double>       m_section_length_vectors;
		std::vector<double>       m_frame_distance_to_section_start;

		// ── Data members (storage only) ───────────────────────────────────
		// The derived class is responsible for calling initData(...) on these
		// from its constructor; CosseratBeamGeometry cannot do it itself
		// because it does not inherit from BaseObject.
		sofa::Data<sofa::type::vector<double>> d_curv_abs_section;
		sofa::Data<sofa::type::vector<double>> d_curv_abs_frames;
		sofa::Data<bool>                       d_debug;

		// ── Internal performance helpers ──────────────────────────────────
		struct JacobianStats {
			double computation_time = 0.0;
			double numerical_error  = 0.0;
			int evaluations_count   = 0;
			int cache_hits          = 0;
			std::chrono::steady_clock::time_point start_time;

			void reset() {
				computation_time = 0.0;
				numerical_error  = 0.0;
				evaluations_count = 0;
				cache_hits = 0;
			}
			void startTiming() { start_time = std::chrono::steady_clock::now(); }
			void endTiming() {
				auto end = std::chrono::steady_clock::now();
				std::chrono::duration<double> elapsed = end - start_time;
				computation_time += elapsed.count();
				evaluations_count++;
			}
		};

		mutable JacobianStats m_jacobian_stats;
		bool m_parallel_computation_enabled = false;

	public:
		CosseratBeamGeometry() = default;
		virtual ~CosseratBeamGeometry() = default;

		// ── Geometry building blocks ──────────────────────────────────────

		/**
		 * @brief Build per-frame topology from `d_curv_abs_section` and
		 *        `d_curv_abs_frames`. Populates `m_frame_to_section_indices*`
		 *        and `m_frame_distance_to_section_start`. Must be called
		 *        every time the discretisation changes.
		 */
		void updateGeometryInfo();

		/**
		 * @brief Build `m_section_properties` from the curvilinear section
		 *        abscissas and the current strain values.
		 */
		void initializeSectionProperties(const StrainVecCoord &strain);

		/**
		 * @brief Compute per-frame arc-length offsets within their section.
		 *        Run after `updateGeometryInfo()`.
		 */
		void initializeFrameProperties();

		/**
		 * @brief Refresh the tangent-exp adjoint matrices on every section
		 *        and frame after strains/transformations have been updated.
		 */
		void updateTangExpSE3();

		// ── Accessors ─────────────────────────────────────────────────────
		const std::vector<SectionInfo> &getSectionProperties() const { return m_section_properties; }
		const std::vector<FrameInfo>   &getFrameProperties()   const { return m_frame_properties; }

		size_t getNumberOfSections() const { return m_section_properties.size(); }
		size_t getNumberOfFrames()   const { return m_frame_properties.size(); }

		void addSection(const SectionInfo &section) { m_section_properties.push_back(section); }
		void addFrame(const FrameInfo &frame)       { m_frame_properties.push_back(frame); }

		void clearSections() { m_section_properties.clear(); }
		void clearFrames()   { m_frame_properties.clear(); }

		void setBeamTopology(const BeamTopology &topology) {
			if (topology.isValid()) {
				m_topology = topology;
			} else {
				msg_warning("CosseratBeamGeometry") << "Invalid beam topology provided";
			}
		}

		const BeamTopology &getBeamTopology() const { return m_topology; }
		bool supportsMultiSectionBeams() const { return true; }

		void enableParallelComputation(bool enable = true) { m_parallel_computation_enabled = enable; }
		bool isParallelComputationEnabled() const { return m_parallel_computation_enabled; }
		void clearComputationCache() { /* placeholder */ }

		// ── Validation utilities ──────────────────────────────────────────
		bool validateSectionProperties() const {
			for (size_t i = 0; i < m_section_properties.size(); ++i) {
				const auto &section = m_section_properties[i];
				if (section.getLength() < 0) {
					msg_warning("CosseratBeamGeometry")
						<< "Section " << i << " has invalid length: " << section.getLength();
					return false;
				}
				if (section.getIndex0() >= section.getIndex1()) {
					msg_warning("CosseratBeamGeometry")
						<< "Section " << i << " has invalid indices: "
						<< section.getIndex0() << " >= " << section.getIndex1();
					return false;
				}
			}
			return true;
		}

		bool checkContinuity(double eps = 1e-6) const {
			for (size_t i = 0; i < m_section_properties.size() - 1; ++i) {
				SE3Type end_transform   = m_section_properties[i].getLocalTransformation(1.0);
				SE3Type start_transform = m_section_properties[i + 1].getLocalTransformation(0.0);
				if (!end_transform.computeIsApprox(start_transform, eps)) {
					msg_warning("CosseratBeamGeometry")
						<< "Discontinuity detected between sections " << i << " and " << i + 1;
					return false;
				}
			}
			return true;
		}

		// ── Trajectory generation ─────────────────────────────────────────
		std::vector<SE3Type> generateSmoothTrajectory(int num_points = 10) const {
			std::vector<SE3Type> trajectory;
			trajectory.reserve(m_section_properties.size() * num_points);
			for (const auto &section : m_section_properties) {
				for (int i = 0; i < num_points; ++i) {
					double t = double(i) / double(num_points);
					trajectory.push_back(section.getLocalTransformation(t));
				}
			}
			if (!m_section_properties.empty()) {
				trajectory.push_back(m_section_properties.back().getLocalTransformation(1.0));
			}
			return trajectory;
		}

		std::vector<SectionInfo> generateSectionTrajectory(int num_points = 10) const {
			std::vector<SectionInfo> trajectory;
			trajectory.reserve(m_section_properties.size() * num_points + 1);
			for (const auto &section : m_section_properties) {
				for (int i = 0; i < num_points; ++i) {
					double t = double(i) / double(num_points);
					SE3Type local_transform = section.getLocalTransformation(t);
					TangentVector current_strain = section.getStrainsVec();
					trajectory.emplace_back(section.getLength(), current_strain, section.getIndex0(), local_transform);
				}
			}
			if (!m_section_properties.empty()) {
				const auto &last_section = m_section_properties.back();
				trajectory.emplace_back(last_section.getLength(), last_section.getStrainsVec(),
										last_section.getIndex0(), last_section.getLocalTransformation(1.0));
			}
			return trajectory;
		}

		// ── Adjoint / tangent matrices ────────────────────────────────────
		std::vector<TangentVector> computeInternalForces(const std::vector<TangentVector> &strains) const {
			std::vector<TangentVector> forces;
			forces.reserve(m_section_properties.size());
			for (size_t i = 0; i < m_section_properties.size(); ++i) {
				if (i < strains.size()) {
					const AdjointMatrix &adj = m_section_properties[i].getAdjoint();
					forces.push_back(adj.transpose() * strains[i]);
				}
			}
			return forces;
		}

		bool validateJacobianAccuracy(double tolerance = 1e-6) const {
			bool all_valid = true;
			for (size_t i = 0; i < m_section_properties.size(); ++i) {
				const auto &section = m_section_properties[i];
				double curv_abs = section.getLength();
				TangentVector strain = section.getStrainsVec();

				AdjointMatrix analytical_jac;
				computeTangExpImplementation(curv_abs, strain, analytical_jac);

				AdjointMatrix numerical_jac = AdjointMatrix::Zero();
				double eps = 1e-7;
				for (int k = 0; k < 6; ++k) {
					TangentVector strain_plus  = strain;
					TangentVector strain_minus = strain;
					strain_plus[k]  += eps;
					strain_minus[k] -= eps;

					SE3Type g_plus  = SE3Type::computeExp(strain_plus * curv_abs);
					SE3Type g_minus = SE3Type::computeExp(strain_minus * curv_abs);
					SE3Type g_diff = g_plus * g_minus.computeInverse();
					TangentVector diff = g_diff.log();
					numerical_jac.col(k) = diff / (2.0 * eps);
				}

				double max_error = (analytical_jac - numerical_jac).cwiseAbs().maxCoeff();
				if (max_error > tolerance) {
					msg_warning("CosseratBeamGeometry")
						<< "Jacobian accuracy check failed for section " << i
						<< ". Max error: " << max_error;
					all_valid = false;
				}
			}
			return all_valid;
		}

		// ── Free-function helpers (static, namespace-level use is fine) ───
		static Matrix3 getTildeMatrix(const Vector3 &u);
		static void buildAdjoint(const Matrix3 &A, const Matrix3 &B, AdjointMatrix &Adjoint);
		static void computeTangExpImplementation(const double &curv_abs, const TangentVector &strain,
												 AdjointMatrix &tang_adjoint_matrix);

	protected:
		// ── Container management ──────────────────────────────────────────
		void reserveContainers(size_t frame_count) {
			m_frame_to_section_indices.clear();
			m_frame_to_section_indices.reserve(frame_count);
			m_frame_to_section_indices_draw.clear();
			m_frame_to_section_indices_draw.reserve(frame_count);
			m_frame_distance_to_section_start.clear();
			m_frame_distance_to_section_start.reserve(frame_count);
		}

		struct SectionIndexResult {
			size_t index_for_frame;
			size_t index_for_next;
			bool   found_exact_match;
		};

		SectionIndexResult findSectionIndex(double frame_curv_abs,
											const sofa::type::vector<double> &sections_curv_abs,
											size_t current_index, double tolerance) {
			if (current_index < sections_curv_abs.size()) {
				const double section_curv_abs = sections_curv_abs[current_index];
				if (std::abs(frame_curv_abs - section_curv_abs) < tolerance) {
					return {current_index, current_index + 1, true};
				} else if (frame_curv_abs > section_curv_abs) {
					return {current_index + 1, current_index + 1, false};
				}
			}
			return {current_index, current_index, false};
		}

		void updateFrameData(size_t frame_idx, size_t section_idx, double frame_curv_absc,
							 const sofa::type::vector<double> &sections_curv_abs) {
			m_frame_to_section_indices.emplace_back(section_idx);
			m_frame_to_section_indices_draw.emplace_back(section_idx);
			m_frame_properties[frame_idx].setRelatedSectionIndex(section_idx);

			const double distance = frame_curv_absc - sections_curv_abs[section_idx - 1];
			m_frame_distance_to_section_start.emplace_back(distance);
			m_frame_properties[frame_idx].setDistanceToSectionStart(std::abs(distance));
		}

		void logCompletionInfo() const {
			if constexpr (ENABLE_GEOMETRY_LOGGING) {
				msg_info("CosseratBeamGeometry")
					<< "updateGeometryInfo completed: m_frame_to_section_indices size = "
					<< m_frame_to_section_indices.size()
					<< ", m_frame_distance_to_section_start size = "
					<< m_frame_distance_to_section_start.size();
			}
		}
	};

} // namespace Cosserat::mapping

#include <Cosserat/mapping/CosseratBeamGeometry.inl>
