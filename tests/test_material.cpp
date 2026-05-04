#include "core/params.hpp"
#include "physics/material.hpp"
#include "physics/soft_rod.hpp"

#include <cassert>
#include <cmath>

int main() {
  EelParams p;
  p.useAspectRatioGeometry = false;
  p.bodyMaterial = BodyMaterial::DragonSkin20;
  p.fluidDensityKgM3 = 1000.0;
  p.physicalBodyLengthM = 0.30;
  p.bodyRadius = 6.0;
  p.eelScale = 150.0;

  const MaterialProperties m = resolveMaterialProperties(p);
  assert(m.name == "dragon_skin_20");
  assert(std::abs(m.densityRatio - 1.08) < 1e-12);
  assert(std::abs(m.modulus100Pa - psiToPa(49.0)) < 1e-6);
  assert(m.youngModulusPa > 5.0e5 && m.youngModulusPa < 6.5e5);
  assert(m.shearModulusPa > 1.8e5 && m.shearModulusPa < 2.1e5);
  assert(m.bulkModulusPa > 9.0e6);
  assert(std::abs(m.tensileStrengthPa - psiToPa(550.0)) < 1e-6);
  assert(std::abs(m.tearStrengthNPerM - pliToNPerM(120.0)) < 1e-9);

  const PlanarRodSectionEstimate rod = estimatePlanarRodSection(p, m);
  assert(rod.valid);
  assert(rod.widthM > 0.0);
  assert(std::abs(rod.thicknessM - 0.02) < 1e-12);
  assert(rod.areaM2 > 0.0);
  assert(rod.secondMomentM4 > 0.0);
  assert(rod.bendingStiffnessNm2 > 0.0);
  assert(rod.axialStiffnessN > 0.0);

  p.rhoBodyRatioOverride = 1.23;
  p.youngModulusPaOverride = 7.5e5;
  p.poissonRatioOverride = 0.45;
  const MaterialProperties overrideMaterial = resolveMaterialProperties(p);
  assert(std::abs(overrideMaterial.densityRatio - 1.23) < 1e-12);
  assert(std::abs(overrideMaterial.youngModulusPa - 7.5e5) < 1e-12);
  assert(std::abs(overrideMaterial.poissonRatio - 0.45) < 1e-12);

  p.physicalBodyLengthM = 0.0;
  const PlanarRodSectionEstimate inactive = estimatePlanarRodSection(p, overrideMaterial);
  assert(!inactive.valid);
}
