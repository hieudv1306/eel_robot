#pragma once

#include "core/params.hpp"
#include "core/types.hpp"
#include "physics/material.hpp"

struct PlanarRodSectionEstimate {
  bool valid = false;
  T physicalLengthM = 0.0;
  T widthM = 0.0;
  T thicknessM = 0.0;
  T areaM2 = 0.0;
  T secondMomentM4 = 0.0;
  T massPerLengthKgM = 0.0;
  T axialStiffnessN = 0.0;
  T bendingStiffnessNm2 = 0.0;
  T dampingRatio = 0.0;
};

PlanarRodSectionEstimate estimatePlanarRodSection(
    const EelParams& p,
    const MaterialProperties& material);
