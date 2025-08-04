#pragma once

#include <Cosserat/config.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>
#include <liegroups/Types.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/core/objectmodel/BaseContext.h>


namespace Cosserat::mapping {

	// Types communs du namespace
	using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
	using SO3Type = sofa::component::cosserat::liegroups::SO3<double>;
	using Vector3 = typename SE3Type::Vector3;
	//using TangentVector = typename SE3Type::TangentVector; // SE3 utilise TangentVector pour les vecteurs 6D
	using Matrix3 = typename SE3Type::Matrix3;
	using Matrix4 = typename SE3Type::Matrix4;
	using AdjointMatrix = typename SE3Type::AdjointMatrix;
	using JacobianMatrix = typename SE3Type::JacobianMatrix;
	using TangentVector = typename SE3Type::TangentVector;

	/**
	 * @brief Class encapsulating the properties of a Cosserat beam node
	 * @todo : change this class to node info instead of section info
	 */
	class SectionInfo {
	private:
		double sec_length_ = 0.0;
		Vector3 angular_strain_ = Vector3::Zero(); // angular strain
		Vector3 linear_strain_ = Vector3::Zero(); // linear strain
		TangentVector strain_ = TangentVector::Zero() ; // strain_ = (angular_strain_^T, linear_strain_^T)^T

		unsigned int index_0_ = 0;
		unsigned int index_1_ = 1;

		// Transformation SE3 au lieu de Matrix4 simple
		SE3Type gX_;

		//Matrix computed automatically
		mutable AdjointMatrix adjoint_;
		mutable AdjointMatrix coAdjoint_;
		mutable bool adjoint_computed_ = false;
		AdjointMatrix tang_adjoint_;
		AdjointMatrix tang_co_adjoint_;

	public:
		SectionInfo() = default;

		// Constructor for node info with length and strain
		// Constructor for node info with length and strain
		// SectionInfo(double length, const Vector3 &kappa, const unsigned int i0,
		// 			const SE3Type &transform = SE3Type::computeIdentity()) :
		// 	sec_length_(length), angular_strain_(kappa), index_0_(i0), gX_(transform) {}

		// Constructor for node info with length and strain
		SectionInfo(double length, const TangentVector &strain, const unsigned int i0,
					const SE3Type &transform = SE3Type::computeIdentity()) :
			sec_length_(length), strain_(strain), index_0_(i0), gX_(transform) {}

		// Constructor for node info with length, strain and index
		SectionInfo(double length, const TangentVector &strain, const unsigned int i0) :
			sec_length_(length), strain_(strain), index_0_(i0) {}

		SectionInfo(double length, const Vector3 &angular_strain, const unsigned int i0) :
			sec_length_(length), angular_strain_(angular_strain), index_0_(i0) {}

		// Accesseurs de base
		double getLength() const { return sec_length_; }
		void setLength(double length) {
			// Verify the length of the beam !!!
			std::cout << "Setting section length to: " << length << std::endl;
			if (length <= 0)
				throw std::invalid_argument("_Length must be positive_");
			sec_length_ = length;
		}
		template<typename VecType>
		void setStrain(const VecType &strain) {
			// Handle SOFA Vec types which use size() instead of SizeAtCompileTime
			if constexpr (std::is_same_v<VecType, sofa::type::Vec<3, double>>) {
				// For sofa::type::Vec<3, double>, convert to our Vector3 type
				for (int i = 0; i < 3; ++i) {
					angular_strain_[i] = strain[i];
				}
				strain_.head<3>() = angular_strain_;
				strain_.tail<3>() = linear_strain_;
			} else if constexpr (std::is_same_v<VecType, sofa::type::Vec<6, double>>) {
				// For sofa::type::Vec<6, double>, split into kappa and gamma
				for (int i = 0; i < 3; ++i) {
					angular_strain_[i] = strain[i];
					linear_strain_[i] = strain[i + 3];

				}
				strain_.head<3>() = angular_strain_;
				strain_.tail<3>() = linear_strain_;
			} else {
				// For Eigen-like types with SizeAtCompileTime
				if constexpr (requires { VecType::SizeAtCompileTime; }) {
					if constexpr (VecType::SizeAtCompileTime == 3) {
						angular_strain_ = strain;
					} else if constexpr (VecType::SizeAtCompileTime == 6) {
						angular_strain_ = strain.template head<3>();
						linear_strain_ = strain.template tail<3>();
					}
					strain_.head<3>() = angular_strain_;
					strain_.tail<3>() = linear_strain_;
				} else {
					// Fallback: use size() method for runtime size determination
					if (strain.size() == 3) {
						for (int i = 0; i < 3; ++i) {
							angular_strain_[i] = strain[i];
						}
					} else if (strain.size() == 6) {
						for (int i = 0; i < 3; ++i) {
							angular_strain_[i] = strain[i];
							linear_strain_[i] = strain[i + 3];
						}
					}
					strain_.head<3>() = angular_strain_;
					strain_.tail<3>() = linear_strain_;
				}
			}
		}

		const TangentVector &getStrainsVec() const { return strain_;}
		// void setStrain(const sofa::type::Vec3 &k) {for (unsigned int i = 0; i<3; i++) kappa_[i] = k[i];}
		//
		// void setStrain(const TangentVector &strain) {
		// 	for (unsigned int i = 0; i < 3; i++) {
		// 		kappa_[i] = strain[i];
		// 		gamma_[i] = strain[i + 3];
		// 	}
		// }
		// void setStrain(const Vector3 &strain) { kappa_ = strain; }

		unsigned int getIndex0() const { return index_0_; }
		unsigned int getIndex1() const { return index_1_; }
		void setIndices(unsigned int i0) {
			index_0_ = i0;
		}

		// Accesseurs pour la transformation SE3
		const SE3Type &getTransformation() const { return gX_; }
		void setTransformation(const SE3Type &transform) {
			gX_ = transform;
			adjoint_computed_ = false; // Invalider le cache
		}

		// Méthodes exploitant les fonctionnalités SE3
		Matrix4 getTransformationMatrix() const { return gX_.matrix(); }

		void setTransformationFromMatrix(const Matrix4 &matrix) {
			gX_ = SE3Type(matrix);
			adjoint_computed_ = true;
		}

		// Exploitation de computeAdjoint() avec cache pour les performances
		const AdjointMatrix &getAdjoint() const {
			if (!adjoint_computed_) {
				adjoint_ = gX_.computeAdjoint();
				coAdjoint_ = adjoint_.transpose();
				adjoint_computed_ = true;
			}
			return adjoint_;
		}

		const AdjointMatrix &getCoAdjoint() const {
			if (!adjoint_computed_) {
				adjoint_ = getAdjoint(); // Compute adjoint and co-adjoint matrix
				coAdjoint_ = adjoint_.transpose();
				return  coAdjoint_;
			}
			return coAdjoint_;
		}

		const AdjointMatrix & getTangAdjointMatrix() {
			return tang_adjoint_;
		}

		void setTanAdjointMatrix(const AdjointMatrix & tang_adjoint_mat) {
				tang_adjoint_ = tang_adjoint_mat;
		}



		// Nouvelles méthodes exploitant les fonctionnalités SE3

		/**
		 * @brief Calcule la transformation locale le long de la section
		 * @param t Paramètre local [0,1] le long de la section
		 * @return Transformation SE3 à la position t
		 */
		SE3Type getLocalTransformation(double t) const {
			if (t < 0.0 || t > 1.0) {
				throw std::invalid_argument("Parameter t must be in [0,1]");
			}

			// Utiliser l'exponentielle SE3 pour calculer la transformation locale
			TangentVector xi;
			xi.template head<3>() = t * sec_length_ * Vector3::UnitX(); // Translation le long de x
			xi.template tail<3>() = t * angular_strain_; // Rotation basée sur la courbure

			return gX_ * SE3Type::computeExp(xi);
		}

		/**
		 * @brief Calcule la dérivée de la transformation par rapport au paramètre t
		 * @param t Paramètre local [0,1]
		 * @return Vecteur tangent représentant la dérivée
		 */
		TangentVector getTransformationDerivative(double /*t*/) const {
			TangentVector xi;
			xi.template head<3>() = sec_length_ * Vector3::UnitX(); // Vitesse de translation
			xi.template tail<3>() = angular_strain_; // Vitesse angulaire

			// Appliquer l'adjoint pour transformer dans le bon repère
			return getAdjoint() * xi;
		}

		/**
		 * @brief Calcule la distance entre deux configurations de section
		 * @param other Autre section à comparer
		 * @param w_rot Poids pour la rotation
		 * @param w_trans Poids pour la translation
		 * @return Distance pondérée
		 */
		double distanceTo(const SectionInfo &other, double w_rot = 1.0, double w_trans = 1.0) const {
			return gX_.distance(other.gX_, w_rot, w_trans);
		}

		/**
		 * @brief Interpolation entre deux sections
		 * @param other Section cible
		 * @param t Paramètre d'interpolation [0,1]
		 * @return Section interpolée
		 * @todo re-implement
		 */
		// SectionInfo interpolate(const SectionInfo &other, double t) const {
		// 	SE3Type interp_transform = gX_.interpolate(other.gX_, t);
		// 	Vector3 interp_kappa = (1.0 - t) * angular_strain_ + t * other.angular_strain_;
		// 	double interp_length = (1.0 - t) * sec_length_ + t * other.sec_length_;
		//
		// 	return SectionInfo(interp_length, interp_kappa, index_0_, interp_transform);
		// }

		/**
		 * @brief Applique une action SE3 à un point
		 * @param point Point à transformer
		 * @return Point transformé
		 */
		Vector3 transformPoint(const Vector3 &point) const { return gX_.computeAction(point); }

		/**
		 * @brief Vérifie si deux sections sont approximativement égales
		 * @param other Autre section
		 * @param eps Tolérance
		 * @return true si les sections sont approximativement égales
		 */
		bool isApprox(const SectionInfo &other, double eps = 1e-6) const {
			return gX_.computeIsApprox(other.gX_, eps) &&
				   (angular_strain_ - other.angular_strain_).norm() < eps &&
				   std::abs(sec_length_ - other.sec_length_) < eps;
		}

		/**
		 * @brief Calcule l'inverse de la transformation
		 * @return Section avec transformation inverse
		 */
		SectionInfo inverse() const {
			return SectionInfo(sec_length_, -strain_, index_0_, gX_.computeInverse());
		}

		/**
		 * @brief Compose deux sections
		 * @param other Section à composer
		 * @return Section composée
		 */
	SectionInfo compose(const SectionInfo &other) const {
		SE3Type composed_transform = gX_.compose(other.gX_);
		// Create a proper 6D strain vector by combining angular and linear strains
		TangentVector composed_strain;
		composed_strain.head<3>() = angular_strain_ + other.angular_strain_; // Angular strain
		composed_strain.tail<3>() = linear_strain_ + other.linear_strain_; // Linear strain
		double total_length = sec_length_ + other.sec_length_;

		return SectionInfo(total_length, composed_strain, index_0_, composed_transform);
	}
	};

	/**
	 * @brief Classe pour les propriétés des frames (utilise TangentVector pour kappa)
	 */
	class FrameInfo {
	private:
		double frames_sect_length_ = 0.0;
		// For frames, we use TangentVector to allow for full curvature representation
		// angular strain (kappa) and linear strain (q)
		TangentVector kappa_ = TangentVector::Zero();
		unsigned int index_0_ = 0;
		unsigned int index_1_ = 1;
		unsigned int related_beam_index_ = 0; // Index de la tige associée
		double distance_to_nearest_beam_node = 0.0; // The distance to the nearest beam node from the base
		SE3Type transformation_;

		mutable AdjointMatrix adjoint_;
		mutable AdjointMatrix coAdjoint_;
		mutable bool adjoint_computed_ = false;
		AdjointMatrix tang_adjoint_;

	public:
		FrameInfo() = default;

		// Utiles methods
		double getLength() const { return frames_sect_length_; }
		void setLength(double length) {
			if (length <= 0)
				throw std::invalid_argument("Length must be positive");
			frames_sect_length_ = length;
		}

		unsigned int get_related_beam_index_() const { return related_beam_index_; }
		void set_related_beam_index_(unsigned int index) {
			if (index < 0)
				throw std::invalid_argument("related_beam_index_ must be non-negative");
			related_beam_index_ = index;
		}

		double getDistanceToNearestBeamNode() const { return distance_to_nearest_beam_node; }
		void setDistanceToNearestBeamNode(double distance) {
			if (distance < 0)
				throw std::invalid_argument("Distance to nearest beam node must be non-negative");
			distance_to_nearest_beam_node = distance;
		}

		const TangentVector &getKappa() const { return kappa_; }
		void setKappa(const TangentVector &k) { kappa_ = k; }

		const SE3Type &getTransformation() const { return transformation_; }
		void setTransformation(const SE3Type &transform) {
			transformation_ = transform;

			//Do I really need this?
			adjoint_computed_ = false;
		}

		const AdjointMatrix &getAdjoint() const {
			if (!adjoint_computed_) {
				adjoint_ = transformation_.computeAdjoint();
				coAdjoint_ = adjoint_.transpose();
				adjoint_computed_ = true;
			}
			return adjoint_;
		}

		const AdjointMatrix & getTangAdjointMatrix() {
			return tang_adjoint_;
		}


		void setTanAdjointMatrix(const AdjointMatrix & tang_adjoint_mat) {
			tang_adjoint_ = tang_adjoint_mat;
		}

		/**
		 * @brief Calcule la transformation locale complète (6D)
		 * @param t Paramètre local [0,1]
		 * @return Transformation SE3
		 */
		SE3Type getLocalTransformation(double t) const {
			if (t < 0.0 || t > 1.0) {
				throw std::invalid_argument("Parameter t must be in [0,1]");
			}

			// Utiliser directement le vecteur kappa_ 6D
			TangentVector xi = t * kappa_;
			xi.template head<3>() += t * frames_sect_length_ * Vector3::UnitX();

			return transformation_ * SE3Type::computeExp(xi);
		}




		/**
		 * @brief Stream output operator for FrameInfo
		 * @param os Output stream
		 * @param frame FrameInfo object to output
		 * @return Reference to output stream
		 */
		friend std::ostream& operator<<(std::ostream& os, const FrameInfo& frame) {
			os << "FrameInfo{length=" << frame.frames_sect_length_
			   << ", related_beam_index=" << frame.related_beam_index_
			   << ", distance_to_nearest=" << frame.distance_to_nearest_beam_node
			   << ", kappa=[" << frame.kappa_.transpose() << "]"
			   << ", transformation=" << frame.transformation_ << "}";
			return os;
		}
	};

	template<class TIn1, class TIn2, class TOut>
	class HookeSeratBaseMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut> {
	public:
		SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE3(HookeSeratBaseMapping, TIn1, TIn2, TOut),
							SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut));


		static constexpr bool ENABLE_GEOMETRY_LOGGING = true;
		using In1 = TIn1;
		using In2 = TIn2;
		using Out = TOut;
		using Inherit = sofa::core::Multi2Mapping<TIn1, TIn2, TOut>;

	protected:
		std::vector<SectionInfo> m_section_properties;
		std::vector<FrameInfo> m_frameProperties;

		// This should be changed by the new Data
		// Geometry information vectors (similar to BaseCosseratMapping)
		std::vector<unsigned int> m_indices_vectors;
		std::vector<unsigned int> m_indices_vectors_draw;
		std::vector<double> m_beam_length_vectors;
		// Additional geometry vectors (like BaseCosseratMapping)
		std::vector<double> m_frames_length_vectors;

		// Helper methods for initialization
		void updateGeometryInfo();
		void initializeSectionProperties();
		void initializeFrameProperties();


		// Validation exploitant les nouvelles méthodes
		bool validateSectionProperties() const {
			for (size_t i = 0; i < m_section_properties.size(); ++i) {
				const auto &section = m_section_properties[i];
				if (section.getLength() < 0) {
					msg_warning() << "Section " << i << " has invalid length: " << section.getLength();
					return false;
				}
				if (section.getIndex0() >= section.getIndex1()) {
					msg_warning() << "Section " << i << " has invalid indices: " << section.getIndex0()
								  << " >= " << section.getIndex1();
					return false;
				}
			}
			return true;
		}

		/**
		 * @brief Calcule la continuité entre sections adjacentes
		 * @param eps Tolérance pour la continuité
		 * @return true si toutes les sections sont continues
		 */
		bool checkContinuity(double eps = 1e-6) const {
			for (size_t i = 0; i < m_section_properties.size() - 1; ++i) {
				SE3Type end_transform = m_section_properties[i].getLocalTransformation(1.0);
				SE3Type start_transform = m_section_properties[i + 1].getLocalTransformation(0.0);

				if (!end_transform.computeIsApprox(start_transform, eps)) {
					msg_warning() << "Discontinuity detected between sections " << i << " and " << i + 1;
					return false;
				}
			}
			return true;
		}

		/**
		 * @brief Génère une trajectoire lisse le long de la tige
		 * @param num_points Nombre de points par section
		 * @return Vecteur de transformations SE3
		 */
		std::vector<SE3Type> generateSmoothTrajectory(int num_points = 10) const {
			std::vector<SE3Type> trajectory;

			for (const auto &section: m_section_properties) {
				for (int i = 0; i < num_points; ++i) {
					double t = double(i) / double(num_points);
					trajectory.push_back(section.getLocalTransformation(t));
				}
			}

			return trajectory;
		}

		/**
		 * @brief Calcule les forces internes utilisant les matrices adjointes
		 * @param strains Déformations d'entrée
		 * @return Forces internes
		 */
		std::vector<TangentVector> computeInternalForces(const std::vector<TangentVector> &strains) const {
			std::vector<TangentVector> forces;
			forces.reserve(m_section_properties.size());

			for (size_t i = 0; i < m_section_properties.size(); ++i) {
				if (i < strains.size()) {
					// Utiliser l'adjoint pour transformer les forces
					const AdjointMatrix &adj = m_section_properties[i].getAdjoint();
					forces.push_back(adj.transpose() * strains[i]);
				}
			}
			return forces;
		}

	public:
		void init() override;
		virtual void doBaseCosseratInit() = 0;
		// void init() override {
		//     try {
		//         // Validation des données
		//         if (!validateSectionProperties()) {
		//             throw std::runtime_error("Invalid section properties");
		//         }
		//
		//         // Vérification de la continuité
		//         if (!checkContinuity()) {
		//             msg_warning() << "Rod sections are not continuous";
		//         }
		//
		//         // Pré-calcul des matrices adjointes (se fait automatiquement avec le cache)
		//         for (const auto& section : m_sectionProperties) {
		//             section.getAdjoint(); // Force le calcul
		//         }
		//
		//         // Initialisation des états
		//         m_strain_state = dynamic_cast<sofa::core::State<In1>*>(this->fromModels1[0].get());
		//         m_rigid_base = dynamic_cast<sofa::core::State<In2>*>(this->fromModels2[0].get());
		//         m_frames = dynamic_cast<sofa::core::State<Out>*>(this->toModels[0].get());
		//
		//         if (!m_strain_state || !m_rigid_base || !m_frames) {
		//             throw std::runtime_error("Failed to initialize mechanical states");
		//         }
		//
		//         Inherit::init();
		//
		//     } catch (const std::exception& e) {
		//         msg_error() << "Initialization failed: " << e.what();
		//         throw;
		//     }
		// }

		// Public methods
		const std::vector<SectionInfo> &getSectionProperties() const { return m_section_properties; }
		const std::vector<FrameInfo> &getFrameProperties() const { return m_frameProperties; }

		size_t getNumberOfSections() const { return m_section_properties.size(); }
		size_t getNumberOfFrames() const { return m_frameProperties.size(); }

		void addSection(const SectionInfo &section) { m_section_properties.push_back(section); }
		void addFrame(const FrameInfo &frame) { m_frameProperties.push_back(frame); }

		void clearSections() { m_section_properties.clear(); }
		void clearFrames() { m_frameProperties.clear(); }

		void updateTangExpSE3();
		//void computeTangExp(double &x, const TangentVector &k, AdjointMatrix &TgX);
		static void computeTangExpImplementation(const double& curv_abs,
	const TangentVector & strain, const AdjointMatrix &adjoint_matrix, AdjointMatrix & tang_adjoint_matrix);

	private:
		struct SectionIndexResult {
			size_t index_for_frame; // Pour set_related_beam_index_
			size_t index_for_next; // Pour la prochaine itération
			bool found_exact_match; // Si on a trouvé une correspondance exacte
		};

		/**
		 * @brief Réserve de l'espace pour les conteneurs de géométrie
		 * @param frame_count Nombre de frames à réserver
		 */
		void reserveContainers(size_t frame_count) {
			m_indices_vectors.clear();
			m_indices_vectors.reserve(frame_count);
			m_indices_vectors_draw.clear();
			m_indices_vectors_draw.reserve(frame_count);
			m_frames_length_vectors.clear();
			m_frames_length_vectors.reserve(frame_count);
		}


		SectionIndexResult findSectionIndex(double frame_curv_abs, const auto &sections_curv_abs, size_t current_index,
											double tolerance) {
			if (current_index < sections_curv_abs.size()) {
				const double section_curv_abs = sections_curv_abs[current_index];

				if (std::abs(frame_curv_abs - section_curv_abs) < tolerance) {
					// Correspondance exacte : utiliser current_index pour la frame, incrémenter pour la suite
					return {current_index, current_index + 1, true};
				} else if (frame_curv_abs > section_curv_abs) {
					return {current_index + 1, current_index + 1, false};
				}
			}
			return {current_index, current_index, false};
		}


		// size_t findSectionIndex(double frame_curv_abs, const auto &sections__curv_abs, size_t current_index,
		//                         double tolerance) {
		//   // Find the section index based on the frame's curvilinear abscissa
		//   if (current_index < sections__curv_abs.size()) {
		//     const double section_curv_abs = sections__curv_abs[current_index];
		//     if (std::abs(frame_curv_abs - section_curv_abs) < tolerance) {
		//       return current_index + 1; // Passage à la section suivante
		//     } else if (frame_curv_abs > section_curv_abs) {
		//       return current_index + 1;
		//     }
		//   }
		//   return current_index;
		// }

		/**
		 * @brief Met à jour les données de frame en function de la section et de la position
		 * @param frame_idx Frame index
		 * @param section_idx section  index
		 * @param frame_curv_absc frame's curvilinear abscissa
		 * @param sections_curv_abs Sections curvilinear abscissas
		 */
		void updateFrameData(size_t frame_idx, size_t section_idx, double frame_curv_absc,
							 const auto &sections_curv_abs) {
			m_indices_vectors.emplace_back(section_idx);
			m_indices_vectors_draw.emplace_back(section_idx);
			m_frameProperties[frame_idx].set_related_beam_index_(section_idx);

			const double distance = frame_curv_absc - sections_curv_abs[section_idx - 1];
			m_frames_length_vectors.emplace_back(distance);
			m_frameProperties[frame_idx].setDistanceToNearestBeamNode(std::abs(distance));
		}
		/**
		 * @brief Log information about the completion of the geometry update
		 */
		// This method is used to log the completion of the geometry update process.
		void logCompletionInfo() const {
			if constexpr (ENABLE_GEOMETRY_LOGGING) { // Constante de compilation
				std::cout<<"HookeSeratBaseMapping updateGeometryInfo completed: m_indices_vectors: " <<
					m_indices_vectors.size() <<std::endl;
				std::cout<< " elements m_frames_length_vectors: " << m_frames_length_vectors.size() << " elements"<<std::endl;
			}
		}

	protected:
		sofa::Data<sofa::type::vector<double>> d_curv_abs_section;
		sofa::Data<sofa::type::vector<double>> d_curv_abs_frames;
		sofa::Data<bool> d_debug;

		// The strain state of the beam, known as \xi in Vec3 or Vec6
		// \xi = (\kappa^T, q^T)^T
		// where \kappa is the angular strain and q is the linear strain
		sofa::core::State<In1> *m_strain_state;
		// The Beam base
		sofa::core::State<In2> *m_rigid_base;
		// The configuration in global frame, known as g(X) in SE(3)
		sofa::core::State<Out> *m_frames;

	protected:
		HookeSeratBaseMapping();
		~HookeSeratBaseMapping() override = default;

		// HookeSeratBaseMapping(const HookeSeratBaseMapping &) = delete;
		HookeSeratBaseMapping &operator=(const HookeSeratBaseMapping &) = delete;
	};

#if !defined(SOFA_COSSERAT_CPP_HookeSeratBaseMapping)
	extern template class SOFA_COSSERAT_API HookeSeratBaseMapping<
			sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
// extern template class SOFA_COSSERAT_API HookeSeratBaseMapping<
// 		sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
