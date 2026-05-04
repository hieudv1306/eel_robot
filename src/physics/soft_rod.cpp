#include "physics/soft_rod.hpp"

#include <cmath>

PlanarRodSectionEstimate estimatePlanarRodSection(
    const EelParams& p,
    const MaterialProperties& material)
{
  PlanarRodSectionEstimate out;
  if (!(p.physicalBodyLengthM > T(0)) ||
      !(p.totalGeometricLengthLU() > T(0)) ||
      !(p.bodyWidthLU() > T(0))) {
    return out;
  }

  out.valid = true;
  out.physicalLengthM = p.physicalBodyLengthM;
  out.widthM = p.bodyWidthLU() / p.totalGeometricLengthLU()
             * p.physicalBodyLengthM;
  out.thicknessM = (p.bodyThicknessM > T(0)) ? p.bodyThicknessM : out.widthM;

  if (!(out.widthM > T(0)) || !(out.thicknessM > T(0))) {
    out.valid = false;
    return out;
  }

  // Rectangular section approximation for a planar soft backbone.  The 2-D
  // capsule width maps to the in-plane bending dimension; thickness is the
  // out-of-plane extrusion used to turn material data into EI.
  out.areaM2 = out.widthM * out.thicknessM;
  out.secondMomentM4 = out.thicknessM * std::pow(out.widthM, T(3)) / T(12);
  out.massPerLengthKgM = material.densityKgM3 * out.areaM2;
  out.axialStiffnessN = material.youngModulusPa * out.areaM2;
  out.bendingStiffnessNm2 = material.youngModulusPa * out.secondMomentM4;
  out.dampingRatio = material.dampingRatio;
  return out;
}
