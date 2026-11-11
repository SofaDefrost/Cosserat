// Micro-benchmark du coût à l'exécution : chaîne de poses SE3 (compose, inverse, log) sur un tableau.
// Mesure le temps par opération et la taille du tableau (le tableau de N poses occupe N*sizeof(SE3)).
#include <chrono>
#include <iostream>
#include <vector>
#include "liegroups/SE3.h"
using namespace sofa::component::cosserat::liegroups;
int main() {
  using G = SE3<double>;
  const int N = 1 << 16, PASSES = 400;
  std::vector<G> g; g.reserve(N);
  for (int i = 0; i < N; ++i) {
    G::TangentVector xi; xi << 0.01*std::sin(i), 0.01*std::cos(i), 0.005*std::sin(2*i), 0.02, 0.01*std::cos(3*i), 0.03;
    g.push_back(G::exp(xi));
  }
  double sink = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int p = 0; p < PASSES; ++p) {
    G acc = G::Identity();
    for (int i = 0; i < N; ++i) { acc = acc * g[i]; if ((i & 1023) == 1023) acc = G::Identity(); }
    sink += acc.translation().sum();
  }
  auto t1 = std::chrono::steady_clock::now();
  for (int p = 0; p < PASSES; ++p)
    for (int i = 0; i < N; ++i) sink += g[i].inverse().translation()[0];
  auto t2 = std::chrono::steady_clock::now();
  for (int p = 0; p < PASSES/4; ++p)
    for (int i = 0; i < N; ++i) sink += g[i].log()[0];
  auto t3 = std::chrono::steady_clock::now();
  auto ns = [](auto a, auto b, double n) { return std::chrono::duration<double, std::nano>(b - a).count() / n; };
  std::cout << "tableau de " << N << " poses : " << (N * sizeof(G)) / 1024 << " Kio (sizeof=" << sizeof(G) << ")\n"
            << "compose : " << ns(t0, t1, double(N) * PASSES) << " ns/op\n"
            << "inverse : " << ns(t1, t2, double(N) * PASSES) << " ns/op\n"
            << "log     : " << ns(t2, t3, double(N) * (PASSES / 4)) << " ns/op   (sink=" << sink << ")\n";
}
