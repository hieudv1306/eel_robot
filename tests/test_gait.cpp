#include "physics/gait.hpp"

#include <cassert>
#include <cmath>

int main() {
  const T eps = 1e-12;
  assert(std::abs(amplitudeEnvelope(0.0, 1.0, 0.07) - 0.007) < eps);
  assert(std::abs(amplitudeEnvelope(1.0, 1.0, 0.07) - 0.161) < eps);

  auto rest = actuationProfile(0.1, 0.2, 1.0);
  assert(rest.tActive == 0.0);
  assert(rest.ampScale == 0.0);
  assert(rest.ampRate == 0.0);

  auto ramp = actuationProfile(0.7, 0.2, 1.0);
  assert(ramp.tActive > 0.0);
  assert(ramp.ampScale > 0.0 && ramp.ampScale < 1.0);
  assert(ramp.ampRate > 0.0);

  auto full = actuationProfile(2.0, 0.2, 1.0);
  assert(full.ampScale == 1.0);
  assert(full.ampRate == 0.0);
}
