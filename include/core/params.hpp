#pragma once

#include "core/enums.hpp"
#include "core/types.hpp"

struct EelParams {
  int    nx           = 2600;
  int    ny           = 400;
  T      initialPositionFactor = 0.85;
  T      inflowVelocity = 0.02;
  T      tau          = 0.55;

  T      bodyRadius   = 6.0;
  T      eelLength    = 1.0;
  T      eelScale     = 150.0;
  int    nSpine       = 200;

  bool   useAspectRatioGeometry = true;
  T      aspectRatio     = 16;
  T      bodyAreaTarget  = 150.0 * 12.0 + 3.14159265358979323846 * 36.0;

  T      eelFreq      = 1.6;
  T      eelLambda    = 1.0;
  T      eelA0        = 0.07;
  GeometryKinematics geometryKinematics = GeometryKinematics::HeightWave;
  T      restTime     = 0.2;
  T      rampTime     = 1.6;

  T      dtAnim       = 0.04;
  int    substeps     = 80;
  T      Ttotal       = 15.0;

  T      kappa        = 2.0;
  int    nIbmIters    = 1;
  int    markerInterval = 1;
  T      warnMeanSlip    = 0.005;
  T      warnMaxSlip     = 0.02;
  T      warnMarkerForce = 0.10;

  T      addedMassFrac = 0.0;

  int    spongeWidth    = 80;
  T      spongeStrength = 0.05;

  int    nWarmup        = 500;

  T omega_lbm() const;
  T nu() const;
  T dtLbm() const;
  int nFrames() const;
  T centerlineLengthLU() const;
  T totalGeometricLengthLU() const;
  T bodyLengthLU() const;
  T bodyWidthLU() const;
  T capsuleAreaLU() const;
  T tailAmplitudeLU() const;
  T tailPeakToPeakLU() const;
};

bool updateGeometryFromAspectRatio(EelParams& p);
