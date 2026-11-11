// Client type cosserat_shell (ShellKinematics.h:261,266,351,424) + SE3::computeAdjoint
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>
using namespace sofa::component::cosserat::liegroups;
int main() {
  Eigen::Vector3d v(0.1, 0.2, 0.3);
  auto a = SO3<double>::computeHat(v);
  auto b = SO3<double>::exp(v).matrix();
  auto c = SO3<double>::computeRightJacobian(v);
  Eigen::Matrix<double, 6, 1> xi; xi << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
  auto g = SE3<double>::exp(0.5 * xi);
  auto d = g.computeAdjoint();
  return (a.norm() + b.norm() + c.norm() + d.norm() > 0) ? 0 : 1;
}
