#include "core/params.hpp"
#include "physics/material.hpp"
#include "physics/soft_backbone.hpp"
#include "physics/soft_rod.hpp"

#include <cassert>
#include <cmath>
#include <numeric>

int main() {
  EelParams p;
  p.useAspectRatioGeometry = false;
  p.bodyMaterial = BodyMaterial::DragonSkin20;
  p.physicalBodyLengthM = 0.30;
  p.bodyThicknessM = 0.02;
  p.bodyRadius = 6.0;
  p.eelScale = 150.0;
  p.nSpine = 80;

  const MaterialProperties material = resolveMaterialProperties(p);
  const PlanarRodSectionEstimate rod = estimatePlanarRodSection(p, material);
  const SoftBackboneConfig config =
    makeSoftBackboneConfig(p, rod, p.nSpine - 1);
  assert(config.valid);
  assert(config.nSegments == p.nSpine - 1);
  assert(config.centerlineLengthM > 0.0);
  assert(std::abs(config.bendingStiffnessNm2 -
                  rod.bendingStiffnessNm2) < 1e-14);
  assert(config.jointAngleStiffnessNm > 0.0);
  assert(config.jointAngleDampingNms > 0.0);

  const SoftBackboneState straight =
    makeStraightBackboneState(config.nSegments);
  SoftBackboneCurvatureProfile zeroPreferred;
  zeroPreferred.curvature.assign(config.nSegments - 1, 0.0);
  zeroPreferred.curvatureRate.assign(config.nSegments - 1, 0.0);
  const SoftBackboneTorqueResult zero =
    computeSoftBackboneTorques(config, straight, zeroPreferred);
  assert(zero.jointMomentNm.size() == static_cast<size_t>(config.nSegments - 1));
  assert(zero.maxAbsJointMomentNm == 0.0);
  assert(zero.elasticEnergyJ == 0.0);

  const T activeTime = p.restTime + p.rampTime + 0.25 / p.eelFreq;
  const SoftBackboneCurvatureProfile preferred =
    preferredBackboneCurvatureWave(p, config, activeTime, p.dtAnim);
  assert(preferred.curvature.size() == static_cast<size_t>(config.nSegments - 1));
  assert(preferred.maxAbsCurvature > 0.0);
  assert(preferred.maxAbsCurvatureRate > 0.0);

  const SoftBackboneState preferredState =
    preferredBackboneStateWave(p, config, activeTime, p.dtAnim);
  assert(preferredState.theta.size() == static_cast<size_t>(config.nSegments));
  assert(preferredState.omega.size() == preferredState.theta.size());

  const SoftBackboneTorqueResult straightDriven =
    computeSoftBackboneTorques(config, straight, preferred);
  assert(straightDriven.maxAbsJointMomentNm > 0.0);
  assert(straightDriven.elasticEnergyJ > 0.0);
  const T torqueSum = std::accumulate(straightDriven.segmentTorqueNm.begin(),
                                      straightDriven.segmentTorqueNm.end(),
                                      T(0));
  assert(std::abs(torqueSum) < 1e-12);

  SoftBackboneState matched =
    makeStraightBackboneState(config.nSegments);
  for (int j = 0; j < config.nSegments - 1; ++j) {
    matched.theta[j + 1] = matched.theta[j] + preferred.curvature[j] * config.dsM;
  }
  SoftBackboneCurvatureProfile preferredNoRate = preferred;
  preferredNoRate.curvatureRate.assign(preferredNoRate.curvature.size(), 0.0);
  const SoftBackboneTorqueResult relaxed =
    computeSoftBackboneTorques(config, matched, preferredNoRate);
  assert(relaxed.maxAbsJointMomentNm < 1e-12);
  assert(relaxed.elasticEnergyJ < 1e-24);

  SoftBackboneState shearing = straight;
  shearing.omega.back() = 1.0;
  const SoftBackboneTorqueResult damped =
    computeSoftBackboneTorques(config, shearing, zeroPreferred);
  assert(damped.maxAbsJointMomentNm > 0.0);
}
