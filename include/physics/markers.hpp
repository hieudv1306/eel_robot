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
