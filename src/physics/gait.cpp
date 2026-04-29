#include "physics/gait.hpp"

#include <algorithm>
#include <cmath>

T amplitudeEnvelope(T s, T L, T A0) {
  T sL = s / L;
  return A0 * (0.1 + 2.2 * sL * sL * sL);
}

ActuationResult actuationProfile(T t, T restTime, T rampTime) {
  T tActive = std::max(T(0), t - restTime);
  if (t < restTime) {
    return {tActive, T(0), T(0)};
  }
  if (rampTime <= T(0)) {
    return {tActive, T(1), T(0)};
  }
  T rampPhase = std::clamp(tActive / rampTime, T(0), T(1));
  T ampScale = rampPhase * rampPhase * (T(3) - T(2) * rampPhase);
  T ampRate = T(0);
  if (rampPhase < T(1)) {
    ampRate = T(6) * rampPhase * (T(1) - rampPhase) / rampTime;
  }
  return {tActive, ampScale, ampRate};
}
