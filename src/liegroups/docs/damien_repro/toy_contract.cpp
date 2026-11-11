// Maquette autonome (Eigen seulement) : deux façons d'imposer le "contrat" CRTP d'un groupe de Lie.
//   (défaut)        CRTP + static_assert    -> aucun coût à l'exécution
//   -DVIRTUAL_FINAL virtual final + hooks purs virtuels (version de Damien)
// Cas d'erreur à compiler un par un :  -DCASE_HIDE  (le dérivé redéfinit inverse())
//                                      -DCASE_MISSING (le dérivé n'implémente pas computeInverse())
#include <Eigen/Geometry>
#include <concepts>
#include <iostream>
#include <type_traits>

#ifdef VIRTUAL_FINAL
// ---- version de Damien : le "private final design pattern" ----
template <class D> struct LieGroup {
  [[nodiscard]] virtual D inverse() const final { return self().computeInverse(); }
  const D& self() const { return static_cast<const D&>(*this); }
 private:
  [[nodiscard]] virtual D computeInverse() const = 0;      // pur virtuel : oblige le dérivé à l'implémenter
};
#else
// ---- alternative : CRTP + vérifications à la compilation, sans virtuel ----
template <class D> struct LieGroup {
  [[nodiscard]] D inverse() const { return self().computeInverse(); }
  const D& self() const { return static_cast<const D&>(*this); }
 protected:
  LieGroup() { contract(); }                               // D est complet ici (corps instancié à l'usage)
 private:
  static constexpr void contract() {
    // 1. hook présent, avec la bonne signature
    static_assert(requires(const D& d) { { d.computeInverse() } -> std::same_as<D>; },
                  "LieGroup: la classe dérivée doit implémenter 'D computeInverse() const'");
    // 2. pas de masquage : &D::inverse doit désigner le membre de la base (type pointeur-membre de LieGroup)
    static_assert(std::is_same_v<decltype(&D::inverse), decltype(&LieGroup::inverse)>,
                  "LieGroup: la classe dérivée ne doit pas redéfinir/masquer inverse()");
  }
};
#endif

struct ToySO3 : LieGroup<ToySO3> {
  Eigen::Quaterniond q{Eigen::Quaterniond::Identity()};
#ifndef CASE_MISSING
  ToySO3 computeInverse() const { ToySO3 r; r.q = q.inverse(); return r; }
#endif
#ifdef CASE_HIDE
  ToySO3 inverse() const { return computeInverse(); }      // le "wrapper" que SO3.h contenait avant 0821e1b
#endif
};
struct ToySE3 : LieGroup<ToySE3> {
  ToySO3 rot; Eigen::Vector3d t{Eigen::Vector3d::Zero()};
  ToySE3 computeInverse() const { ToySE3 r; r.rot = rot.inverse(); r.t = -(r.rot.q * t); return r; }
};

int main() {
  ToySO3 a; ToySE3 b; (void)a.inverse(); (void)b.inverse();
  std::cout << "sizeof(ToySO3)=" << sizeof(ToySO3) << " sizeof(ToySE3)=" << sizeof(ToySE3)
            << " polymorphe=" << std::is_polymorphic_v<ToySE3> << "\n";
}
