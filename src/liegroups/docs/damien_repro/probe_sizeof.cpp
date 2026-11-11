// Taille et polymorphisme des groupes de Lie : à compiler contre l'arbre à tester (-I<arbre>/src).
#include <iostream>
#include <type_traits>
#include "liegroups/SO3.h"
#include "liegroups/SE3.h"
#include "liegroups/SE2.h"
using namespace sofa::component::cosserat::liegroups;
template <class T> void line(const char* n) {
  std::cout << n << " sizeof=" << sizeof(T) << " alignof=" << alignof(T)
            << " polymorphe=" << std::is_polymorphic_v<T> << "\n";
}
int main() {
  line<SO3<double>>("SO3<double>");
  line<SE3<double>>("SE3<double>");
  line<SE2<double>>("SE2<double>");
  line<SO3<float>>("SO3<float> ");
  line<SE3<float>>("SE3<float> ");
  std::cout << "(référence : Eigen::Quaterniond=" << sizeof(Eigen::Quaterniond)
            << ", Vector3d=" << sizeof(Eigen::Vector3d) << ")\n";
}
