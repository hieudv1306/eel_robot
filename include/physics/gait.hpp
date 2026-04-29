#pragma once

#include "core/types.hpp"

struct ActuationResult {
  T tActive;
  T ampScale;
  T ampRate;
};

T amplitudeEnvelope(T s, T L, T A0);
ActuationResult actuationProfile(T t, T restTime, T rampTime);
