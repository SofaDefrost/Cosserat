// Pourquoi (g * g^-1).isApprox(Identity) peut renvoyer faux : quatre manières de comparer à l'identité.
#include <cmath>
#include <iomanip>
#include <iostream>
#include "liegroups/SE3.h"
using namespace sofa::component::cosserat::liegroups;
int main() {
  using G = SE3<double>;
  G::TangentVector xi; xi << 0.1, 0.2, 0.3, 1.0, 2.0, 3.0;     // [phi ; rho]
  const G g = G::exp(xi);
  const G e = g * g.inverse();
  const G I = G::Identity();
  std::cout << std::scientific << std::setprecision(3);
  // 1. API du dépôt : isApprox du groupe (tolérance par défaut ou non selon la version)
  bool api = e.isApprox(I);
  // 2. Ce que fait SE3::computeIsApprox : isApprox d'Eigen composante par composante (test RELATIF)
  const Eigen::Vector3d t = e.translation();
  const double rel_rot = e.rotation().isApprox(I.rotation(), 1e-12) ? 1 : 0;
  const bool rel_t = t.isApprox(Eigen::Vector3d::Zero(), 1e-12);
  // 3. Test absolu : norme de l'écart
  const double abs_t = t.norm();
  const double abs_log = e.log().norm();
  std::cout << "isApprox du groupe (API)                    : " << api << "\n"
            << "rotation.isApprox(I, 1e-12)  (relatif)      : " << rel_rot << "\n"
            << "translation.isApprox(0, 1e-12) (relatif)    : " << rel_t << "   <-- faux si t != 0 exactement\n"
            << "|translation(g*g^-1)|  (absolu)             : " << abs_t << "\n"
            << "|log(g*g^-1)|          (absolu)             : " << abs_log << "\n";
  // Règle d'Eigen : a.isApprox(b, p)  <=>  |a-b|^2 <= p^2 * min(|a|^2, |b|^2)
  // Avec b = 0 le second membre vaut 0 : il faut a == 0 exactement.
  return 0;
}
