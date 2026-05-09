#pragma once

#include "core/params.hpp"
#include "core/types.hpp"
#include "physics/soft_backbone.hpp"

#include <vector>

struct LagrangianMarkers {
  std::vector<T> x, y;
  std::vector<T> ud, vd;
  std::vector<T> ds;
  std::vector<T> fx, fy;

  int size() const;
  void resize(int n);
};

struct SoftBackboneForceProjection {
  std::vector<T> segmentTorqueNm;
  T netForceXLat = 0.0;
  T netForceYLat = 0.0;
  T netTorqueLat = 0.0;
  T maxAbsSegmentTorqueNm = 0.0;
  T latticeTorqueToNm = 0.0;
};

void buildLagrangianMarkers(
    const EelParams& p,
    T t, T Vx, T Vy, T omegaZ,
    T xCm, T yCm, T theta,
    T dtLbm, T ampRamp,
    LagrangianMarkers& markers);

void buildLagrangianMarkersFromSoftBackbone(
    const EelParams& p,
    const SoftBackboneConfig& config,
    T t, T Vx, T Vy, T omegaZ,
    T xCm, T yCm, T theta,
    T dtLbm, T ampRamp,
    LagrangianMarkers& markers);

void buildLagrangianMarkersFromSoftBackboneState(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    T Vx, T Vy, T omegaZ,
    T xCm, T yCm, T theta,
    T dtLbm,
    LagrangianMarkers& markers);

SoftBackboneForceProjection projectMarkerForcesToSoftBackbone(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    T xCm, T yCm, T theta,
    const LagrangianMarkers& markers);
