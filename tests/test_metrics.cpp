#include "physics/diagnostics.hpp"

#include <cassert>
#include <cmath>
#include <vector>

int main() {
  std::vector<T> t{0.0, 0.1, 0.2};
  std::vector<T> u{1.0, 2.0, 3.0};
  std::vector<T> p{2.0, -4.0, 6.0};
  std::vector<T> f{0.5, 0.5, 0.5};
  std::vector<T> lateral{0.0, 0.0, 0.0};
  std::vector<T> re{10.0, 20.0, 30.0};
  std::vector<T> st{0.1, 0.2, 0.3};
  std::vector<T> ustar{0.4, 0.5, 0.6};
  std::vector<T> eta{0.0, 0.0, 0.0};
  std::vector<T> slip{0.01, 0.02, 0.03};
  std::vector<T> maxSlip{0.02, 0.03, 0.04};
  std::vector<T> markerForce{1.0, 1.0, 1.0};
  std::vector<T> maxMarkerForce{2.0, 2.0, 2.0};
  std::vector<T> normSlip{0.1, 0.2, 0.3};
  std::vector<T> powerRigid{1.0, 2.0, 3.0};
  std::vector<T> powerDef{1.0, -6.0, 3.0};
  std::vector<T> residual{0.001, 0.002, 0.003};
  std::vector<T> maxResidual{0.004, 0.005, 0.006};

  auto s = computeSteadySummary(t, u, p, f, lateral, re, st, ustar, eta,
                                slip, maxSlip, markerForce, maxMarkerForce,
                                normSlip, powerRigid, powerDef,
                                residual, maxResidual, 0.0, 2.0);

  assert(s.nSamples == 2);
  assert(std::abs(s.meanU - 2.5) < 1e-12);
  assert(std::abs(s.meanP - 5.0) < 1e-12);
  assert(std::abs(s.meanUstar - 0.55) < 1e-12);
  assert(std::abs(s.meanNormalizedSlip - 0.25) < 1e-12);
  assert(std::abs(s.meanAbsDeformationPower - 4.5) < 1e-12);
  assert(std::abs(s.meanResidualSlip - 0.0025) < 1e-12);
  assert(std::abs(s.maxResidualSlip - 0.006) < 1e-12);
}
