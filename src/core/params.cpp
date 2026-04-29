#include "core/params.hpp"

#include <cmath>

T EelParams::omega_lbm() const { return 1.0 / tau; }
T EelParams::nu() const { return (tau - 0.5) / 3.0; }
T EelParams::dtLbm() const { return dtAnim / substeps; }
int EelParams::nFrames() const { return static_cast<int>(Ttotal / dtAnim); }
T EelParams::centerlineLengthLU() const { return eelScale * eelLength; }
T EelParams::totalGeometricLengthLU() const { return centerlineLengthLU() + 2.0 * bodyRadius; }
T EelParams::bodyLengthLU() const { return centerlineLengthLU(); }
T EelParams::bodyWidthLU() const { return 2.0 * bodyRadius; }
T EelParams::capsuleAreaLU() const {
  return centerlineLengthLU() * bodyWidthLU() + M_PI * bodyRadius * bodyRadius;
}

T EelParams::tailAmplitudeLU() const {
  T env = eelA0 * (0.1 + 2.2);
  return env * eelScale;
}

T EelParams::tailPeakToPeakLU() const { return 2.0 * tailAmplitudeLU(); }

bool updateGeometryFromAspectRatio(EelParams& p)
{
  if (!p.useAspectRatioGeometry) {
    return true;
  }
  if (!(p.aspectRatio > T(1)) || !(p.bodyAreaTarget > T(0))) {
    return false;
  }

  const T denom = p.aspectRatio - T(1) + M_PI / T(4);
  if (!(denom > T(0))) {
    return false;
  }

  const T D  = std::sqrt(p.bodyAreaTarget / denom);
  const T Lc = (p.aspectRatio - T(1)) * D;

  p.bodyRadius = T(0.5) * D;
  p.eelLength  = T(1);
  p.eelScale   = Lc;
  return std::isfinite(p.bodyRadius) && std::isfinite(p.eelScale)
      && p.bodyRadius > T(0) && p.eelScale > T(0);
}
