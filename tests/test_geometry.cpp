#include "core/params.hpp"
#include "physics/geometry.hpp"

#include <algorithm>
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

  std::vector<T> xp, yp, uxp, uyp, dsp;
  std::vector<T> xm, ym, uxm, uym, dsm;
  const T t = 1.0;
  const T dt = p.dtLbm();
  buildCapsuleGeometry(p, t, 100.0, 60.0, 0.0, dt, 1.0,
                       x, y, ux, uy, ds);
  buildCapsuleGeometry(p, t + dt, 100.0, 60.0, 0.0, dt, 1.0,
                       xp, yp, uxp, uyp, dsp);
  buildCapsuleGeometry(p, t - dt, 100.0, 60.0, 0.0, dt, 1.0,
                       xm, ym, uxm, uym, dsm);
  assert(x.size() == xp.size());
  assert(x.size() == xm.size());
  T maxAbsUx = 0.0;
  T maxAbsUy = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    const T uxFd = 0.5 * (xp[i] - xm[i]);
    const T uyFd = 0.5 * (yp[i] - ym[i]);
    assert(std::abs(ux[i] - uxFd) < 1e-12);
    assert(std::abs(uy[i] - uyFd) < 1e-12);
    maxAbsUx = std::max(maxAbsUx, std::abs(ux[i]));
    maxAbsUy = std::max(maxAbsUy, std::abs(uy[i]));
  }
  assert(maxAbsUx > 1e-8);
  assert(maxAbsUy > 1e-8);

  EelParams pInext = p;
  pInext.nSpine = 50;
  pInext.geometryKinematics = GeometryKinematics::InextensibleWave;
  buildCapsuleGeometry(pInext, 2.0, 100.0, 60.0, 0.0, dt, 1.0,
                       x, y, ux, uy, ds);
  const T dsSpine = pInext.eelScale * pInext.eelLength / (pInext.nSpine - 1);
  const int nTailCap = std::max(
    static_cast<int>(std::ceil(3.14159265358979323846 * pInext.bodyRadius / dsSpine)),
    8);
  const int lowerStart = pInext.nSpine + (nTailCap - 1);
  std::vector<T> cx(pInext.nSpine), cy(pInext.nSpine);
  for (int i = 0; i < pInext.nSpine; ++i) {
    const int lowerIdx = lowerStart + (pInext.nSpine - 1 - i);
    cx[i] = T(0.5) * (x[i] + x[lowerIdx]);
    cy[i] = T(0.5) * (y[i] + y[lowerIdx]);
  }
  T centerlineArc = 0.0;
  for (int i = 1; i < pInext.nSpine; ++i) {
    const T dx = cx[i] - cx[i - 1];
    const T dy = cy[i] - cy[i - 1];
    centerlineArc += std::sqrt(dx * dx + dy * dy);
  }
  assert(std::abs(centerlineArc - pInext.centerlineLengthLU()) < 1e-10);

  auto e = bodyForwardAxis(0.3);
  auto n = bodyLateralAxis(0.3);
  const T dot = e[0] * n[0] + e[1] * n[1];
  const T eNorm = std::sqrt(e[0] * e[0] + e[1] * e[1]);
  const T nNorm = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  assert(std::abs(dot) < 1e-12);
  assert(std::abs(eNorm - 1.0) < 1e-12);
  assert(std::abs(nNorm - 1.0) < 1e-12);
}
