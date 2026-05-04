#pragma once

#include "core/params.hpp"
#include "core/types.hpp"

#include <string>

struct MaterialProperties {
  BodyMaterial material = BodyMaterial::DragonSkin20;
  std::string name;
  T shoreA = 0.0;
  T densityKgM3 = 0.0;
  T densityRatio = 1.0;
  T modulus100Pa = 0.0;
  T youngModulusPa = 0.0;
  T poissonRatio = 0.0;
  T shearModulusPa = 0.0;
  T bulkModulusPa = 0.0;
  T tensileStrengthPa = 0.0;
  T elongationAtBreak = 0.0;
  T tearStrengthNPerM = 0.0;
  T dampingRatio = 0.0;
};

T psiToPa(T psi);
T pliToNPerM(T pli);

MaterialProperties defaultDragonSkin20Material(T fluidDensityKgM3);
MaterialProperties neutralLatticeMaterial(T fluidDensityKgM3);
MaterialProperties resolveMaterialProperties(const EelParams& p);
