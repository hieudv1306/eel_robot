#include "core/params.hpp"
#include "physics/geometry.hpp"

#include <cassert>
#include <cmath>
#include <vector>

int main() {
  EelParams p;
  p.nSpine = 20;
  p.bodyRadius = 3.0;
  p.eelScale = 40.0;
  p.useAspectRatioGeometry = false;

  std::vector<T> x, y, ux, uy, ds;
  buildCapsuleGeometry(p, 0.25, 100.0, 60.0, 0.1, p.dtLbm(), 1.0,
                       x, y, ux, uy, ds);

  assert(!x.empty());
  assert(x.size() == y.size());
  assert(x.size() == ux.size());
  assert(x.size() == uy.size());
  assert(x.size() == ds.size());

  auto e = bodyForwardAxis(0.3);
  auto n = bodyLateralAxis(0.3);
  const T dot = e[0] * n[0] + e[1] * n[1];
  const T eNorm = std::sqrt(e[0] * e[0] + e[1] * e[1]);
  const T nNorm = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  assert(std::abs(dot) < 1e-12);
  assert(std::abs(eNorm - 1.0) < 1e-12);
  assert(std::abs(nNorm - 1.0) < 1e-12);
}
