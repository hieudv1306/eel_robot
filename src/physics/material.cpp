#include "physics/material.hpp"

#include <algorithm>
#include <cmath>

namespace {

T clampPoissonRatio(T nu)
{
  return std::clamp(nu, T(-0.99), T(0.499));
}

void fillDerivedElasticConstants(MaterialProperties& m)
{
  m.poissonRatio = clampPoissonRatio(m.poissonRatio);
  if (m.youngModulusPa > T(0)) {
    m.shearModulusPa = m.youngModulusPa / (T(2) * (T(1) + m.poissonRatio));
    m.bulkModulusPa = m.youngModulusPa / (T(3) * (T(1) - T(2) * m.poissonRatio));
  } else {
    m.shearModulusPa = T(0);
    m.bulkModulusPa = T(0);
  }
}

T validFluidDensity(T fluidDensityKgM3)
{
  return (fluidDensityKgM3 > T(0)) ? fluidDensityKgM3 : T(1000);
}

}  // namespace

T psiToPa(T psi)
{
  return psi * T(6894.757293168);
}

T pliToNPerM(T pli)
{
  return pli * T(4.4482216152605) / T(0.0254);
}

MaterialProperties defaultDragonSkin20Material(T fluidDensityKgM3)
{
  MaterialProperties m;
  m.material = BodyMaterial::DragonSkin20;
  m.name = bodyMaterialName(m.material);
  m.shoreA = T(20);
  m.densityKgM3 = T(1080);
  m.densityRatio = m.densityKgM3 / validFluidDensity(fluidDensityKgM3);
  m.modulus100Pa = psiToPa(T(49));

  // Smooth-On publishes the 100% modulus, not a small-strain Young's modulus.
  // For the first soft-rod estimate, treat Dragon Skin 20 as incompressible
  // Neo-Hookean and infer mu from nominal uniaxial stress at lambda = 2:
  //   P = mu * (lambda - lambda^-2).
  const T lambda = T(2);
  const T mu = m.modulus100Pa / (lambda - T(1) / (lambda * lambda));
  m.poissonRatio = T(0.49);
  m.youngModulusPa = T(2) * mu * (T(1) + m.poissonRatio);

  m.tensileStrengthPa = psiToPa(T(550));
  m.elongationAtBreak = T(6.20);
  m.tearStrengthNPerM = pliToNPerM(T(120));
  m.dampingRatio = T(0.05);
  fillDerivedElasticConstants(m);
  return m;
}

MaterialProperties neutralLatticeMaterial(T fluidDensityKgM3)
{
  MaterialProperties m;
  m.material = BodyMaterial::NeutralLattice;
  m.name = bodyMaterialName(m.material);
  m.densityKgM3 = validFluidDensity(fluidDensityKgM3);
  m.densityRatio = T(1);
  m.dampingRatio = T(0);
  return m;
}

MaterialProperties resolveMaterialProperties(const EelParams& p)
{
  MaterialProperties m =
    (p.bodyMaterial == BodyMaterial::NeutralLattice)
      ? neutralLatticeMaterial(p.fluidDensityKgM3)
      : defaultDragonSkin20Material(p.fluidDensityKgM3);

  if (p.rhoBodyRatioOverride > T(0)) {
    m.densityRatio = p.rhoBodyRatioOverride;
    m.densityKgM3 = m.densityRatio * validFluidDensity(p.fluidDensityKgM3);
  }
  if (p.youngModulusPaOverride > T(0)) {
    m.youngModulusPa = p.youngModulusPaOverride;
  }
  if (p.poissonRatioOverride > T(-0.99) && p.poissonRatioOverride < T(0.5)) {
    m.poissonRatio = p.poissonRatioOverride;
  }
  if (p.materialDampingRatio >= T(0)) {
    m.dampingRatio = p.materialDampingRatio;
  }
  fillDerivedElasticConstants(m);
  return m;
}
