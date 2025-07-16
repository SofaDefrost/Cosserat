#pragma once
#include <Cosserat/config.h>

#include <liegroups/SE3.h>
#include <liegroups/SO3.h>
#include <liegroups/Types.h>
#include <sofa/core/Multi2Mapping.h>

namespace Cosserat::mapping {

	// Types communs du namespace
	using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
	using SO3Type = sofa::component::cosserat::liegroups::SO3<double>;
	using Vector3 = typename SE3Type::Vector3;
	using Vector6 = typename SE3Type::TangentVector; // SE3 utilise TangentVector pour les vecteurs 6D
	using Matrix3 = typename SE3Type::Matrix3;
	using Matrix4 = typename SE3Type::Matrix4;
	using AdjointMatrix = typename SE3Type::AdjointMatrix;
	using JacobianMatrix = typename SE3Type::JacobianMatrix;
	using TangentVector = typename SE3Type::TangentVector;

	/**
	 * @brief Classe encapsulant les propriétés d'une section de tige de Cosserat
	 */
	class SectionInfo {
	private:
		double sec_length_ = 0.0;
		Vector3 kappa_ = Vector3::Zero(); // Courbure pour les sections (3D)
		unsigned int index_0_ = 0;
		unsigned int index_1_ = 1;

		// Transformation SE3 au lieu de Matrix4 simple
		SE3Type transformation_;

		// Matrices calculées automatiquement
		mutable AdjointMatrix adjoint_;
		mutable AdjointMatrix coAdjoint_;
		mutable bool adjoint_computed_ = false;

	public:
		SectionInfo() = default;

		// Constructeur avec transformation SE3
		SectionInfo(double length, const Vector3 &kappa, unsigned int i0, unsigned int i1,
					const SE3Type &transform = SE3Type::computeIdentity()) :
			sec_length_(length), kappa_(kappa), index_0_(i0), index_1_(i1), transformation_(transform) {}

		// Accesseurs de base
		double getLength() const { return sec_length_; }
		void setLength(double length) {
			if (length <= 0)
				throw std::invalid_argument("Length must be positive");
			sec_length_ = length;
		}

		const Vector3 &getKappa() const { return kappa_; }
		void setKappa(const Vector3 &k) { kappa_ = k; }

		unsigned int getIndex0() const { return index_0_; }
		unsigned int getIndex1() const { return index_1_; }
		void setIndices(unsigned int i0, unsigned int i1) {
			if (i0 >= i1)
				throw std::invalid_argument("index_0 must be less than index_1");
			index_0_ = i0;
			index_1_ = i1;
		}

		// Accesseurs pour la transformation SE3
		const SE3Type &getTransformation() const { return transformation_; }
		void setTransformation(const SE3Type &transform) {
			transformation_ = transform;
			adjoint_computed_ = false; // Invalider le cache
		}

		// Méthodes exploitant les fonctionnalités SE3
		Matrix4 getTransformationMatrix() const { return transformation_.matrix(); }

		void setTransformationFromMatrix(const Matrix4 &matrix) {
			transformation_ = SE3Type(matrix);
			adjoint_computed_ = false;
		}

		// Exploitation de computeAdjoint() avec cache pour les performances
		const AdjointMatrix &getAdjoint() const {
			if (!adjoint_computed_) {
				adjoint_ = transformation_.computeAdjoint();
				coAdjoint_ = adjoint_.transpose();
				adjoint_computed_ = true;
			}
			return adjoint_;
		}

		const AdjointMatrix &getCoAdjoint() const {
			if (!adjoint_computed_) {
				getAdjoint(); // Compute adjoint and co-adjoint matrix

			}
			return coAdjoint_;
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
			xi.template tail<3>() = t * kappa_; // Rotation basée sur la courbure

			return transformation_ * SE3Type::computeExp(xi);
		}

		/**
		 * @brief Calcule la dérivée de la transformation par rapport au paramètre t
		 * @param t Paramètre local [0,1]
		 * @return Vecteur tangent représentant la dérivée
		 */
		TangentVector getTransformationDerivative(double /*t*/) const {
			TangentVector xi;
			xi.template head<3>() = sec_length_ * Vector3::UnitX(); // Vitesse de translation
			xi.template tail<3>() = kappa_; // Vitesse angulaire

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
			return transformation_.distance(other.transformation_, w_rot, w_trans);
		}

		/**
		 * @brief Interpolation entre deux sections
		 * @param other Section cible
		 * @param t Paramètre d'interpolation [0,1]
		 * @return Section interpolée
		 */
		SectionInfo interpolate(const SectionInfo &other, double t) const {
			SE3Type interp_transform = transformation_.interpolate(other.transformation_, t);
			Vector3 interp_kappa = (1.0 - t) * kappa_ + t * other.kappa_;
			double interp_length = (1.0 - t) * sec_length_ + t * other.sec_length_;

			return SectionInfo(interp_length, interp_kappa, index_0_, index_1_, interp_transform);
		}

		/**
		 * @brief Applique une action SE3 à un point
		 * @param point Point à transformer
		 * @return Point transformé
		 */
		Vector3 transformPoint(const Vector3 &point) const { return transformation_.computeAction(point); }

		/**
		 * @brief Vérifie si deux sections sont approximativement égales
		 * @param other Autre section
		 * @param eps Tolérance
		 * @return true si les sections sont approximativement égales
		 */
		bool isApprox(const SectionInfo &other, double eps = 1e-6) const {
			return transformation_.computeIsApprox(other.transformation_, eps) &&
				   (kappa_ - other.kappa_).norm() < eps && std::abs(sec_length_ - other.sec_length_) < eps;
		}

		/**
		 * @brief Calcule l'inverse de la transformation
		 * @return Section avec transformation inverse
		 */
		SectionInfo inverse() const {
			return SectionInfo(sec_length_, -kappa_, index_1_, index_0_, transformation_.computeInverse());
		}

		/**
		 * @brief Compose deux sections
		 * @param other Section à composer
		 * @return Section composée
		 */
		SectionInfo compose(const SectionInfo &other) const {
			SE3Type composed_transform = transformation_.compose(other.transformation_);
			Vector3 composed_kappa = kappa_ + other.kappa_; // Approximation linéaire
			double total_length = sec_length_ + other.sec_length_;

			return SectionInfo(total_length, composed_kappa, index_0_, other.index_1_, composed_transform);
		}
	};

	/**
	 * @brief Classe pour les propriétés des frames (utilise Vector6 pour kappa)
	 */
	class FrameInfo {
	private:
		double frames_sect_length_ = 0.0;
		Vector6 kappa_ = Vector6::Zero(); // Courbure complète 6D pour les frames
		unsigned int index_0_ = 0;
		unsigned int index_1_ = 1;
		unsigned int beam_index_ = 0; // Index de la tige associée
		SE3Type transformation_;

		mutable AdjointMatrix adjoint_;
		mutable AdjointMatrix coAdjoint_;
		mutable bool adjoint_computed_ = false;

	public:
		FrameInfo() = default;

		// Méthodes similaires à SectionInfo mais adaptées aux frames
		double getLength() const { return frames_sect_length_; }
		void setLength(double length) {
			if (length <= 0)
				throw std::invalid_argument("Length must be positive");
			frames_sect_length_ = length;
		}

		const Vector6 &getKappa() const { return kappa_; }
		void setKappa(const Vector6 &k) { kappa_ = k; }

		const SE3Type &getTransformation() const { return transformation_; }
		void setTransformation(const SE3Type &transform) {
			transformation_ = transform;
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

		// Autres méthodes similaires à SectionInfo...
	};

	template<class TIn1, class TIn2, class TOut>
	class HookeSeratBaseMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut> {
	public:
		SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE3(HookeSeratBaseMapping, TIn1, TIn2, TOut),
							SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut));

		using In1 = TIn1;
		using In2 = TIn2;
		using Out = TOut;
		using Inherit = sofa::core::Multi2Mapping<TIn1, TIn2, TOut>;

	protected:
		std::vector<SectionInfo> m_sectionProperties;
		std::vector<FrameInfo> m_frameProperties;

		// Geometry information vectors (similar to BaseCosseratMapping)
		std::vector<unsigned int> m_indicesVectors;
		std::vector<unsigned int> m_indicesVectorsDraw;
		std::vector<double> m_beamLengthVectors;

		// Helper methods for initialization
		void updateGeometryInfo();
		void initializeSectionProperties();
		void initializeFrameProperties();

		// Validation exploitant les nouvelles méthodes
		bool validateSectionProperties() const {
			for (size_t i = 0; i < m_sectionProperties.size(); ++i) {
				const auto &section = m_sectionProperties[i];
				if (section.getLength() <= 0) {
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
			for (size_t i = 0; i < m_sectionProperties.size() - 1; ++i) {
				SE3Type end_transform = m_sectionProperties[i].getLocalTransformation(1.0);
				SE3Type start_transform = m_sectionProperties[i + 1].getLocalTransformation(0.0);

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

			for (const auto &section: m_sectionProperties) {
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
		std::vector<Vector6> computeInternalForces(const std::vector<Vector6> &strains) const {
			std::vector<Vector6> forces;
			forces.reserve(m_sectionProperties.size());

			for (size_t i = 0; i < m_sectionProperties.size(); ++i) {
				if (i < strains.size()) {
					// Utiliser l'adjoint pour transformer les forces
					const AdjointMatrix &adj = m_sectionProperties[i].getAdjoint();
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
		const std::vector<SectionInfo> &getSectionProperties() const { return m_sectionProperties; }
		const std::vector<FrameInfo> &getFrameProperties() const { return m_frameProperties; }

		size_t getNumberOfSections() const { return m_sectionProperties.size(); }
		size_t getNumberOfFrames() const { return m_frameProperties.size(); }

		void addSection(const SectionInfo &section) { m_sectionProperties.push_back(section); }
		void addFrame(const FrameInfo &frame) { m_frameProperties.push_back(frame); }

		void clearSections() { m_sectionProperties.clear(); }
		void clearFrames() { m_frameProperties.clear(); }

	protected:
    sofa::Data<sofa::type::vector<double>> d_curv_abs_section;
    sofa::Data<sofa::type::vector<double>> d_curv_abs_frames;
		sofa::Data<bool> d_debug;

		sofa::core::State<In1> *m_strain_state;
		sofa::core::State<In2> *m_rigid_base;
		sofa::core::State<Out> *m_frames;

	protected:
		HookeSeratBaseMapping();
		~HookeSeratBaseMapping() override = default;

		HookeSeratBaseMapping(const HookeSeratBaseMapping &) = delete;
		HookeSeratBaseMapping &operator=(const HookeSeratBaseMapping &) = delete;
	};

#if !defined(SOFA_COSSERAT_CPP_HookeSeratBaseMapping)
	extern template class SOFA_COSSERAT_API HookeSeratBaseMapping<
			sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
	extern template class SOFA_COSSERAT_API HookeSeratBaseMapping<
			sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
