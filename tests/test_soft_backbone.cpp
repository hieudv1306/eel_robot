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
  // Added mass: full theoretical fraction (default).  For Dragon Skin 20 in
  // water with the Issue 1 validation params, the slender-body added inertia
  // should be of comparable magnitude to the structural inertia (within an
  // order of magnitude); a zero or negative value would indicate the
  // partitioned-FSI added-mass mitigation has been disabled by accident.
  assert(config.addedSegmentRotationalInertiaKgM2 > 0.0);
  assert(config.effectiveSegmentRotationalInertiaKgM2 >
         config.segmentRotationalInertiaKgM2);
  // Disabling the added-mass fraction must zero out the augmentation.
  const SoftBackboneConfig configStructOnly =
    makeSoftBackboneConfig(p, rod, p.nSpine - 1, T(0));
  assert(configStructOnly.addedSegmentRotationalInertiaKgM2 == 0.0);
  assert(std::abs(configStructOnly.effectiveSegmentRotationalInertiaKgM2 -
                  configStructOnly.segmentRotationalInertiaKgM2) < 1e-30);

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

  // ----- Implicit (2nd-order Newton-Euler) integrator invariants -----
  // (a) Static equilibrium: preferred at rest (zero theta, zero omega) is a
  //     fixed point of the implicit Euler operator with no fluid load.  A
  //     state initially at this rest state must remain at it bit-for-bit.
  SoftBackboneState restingPreferredState =
    makeStraightBackboneState(config.nSegments);
  SoftBackboneState implicitState = restingPreferredState;
  SoftBackboneDynamicsParams implicitParams;
  implicitParams.maxAngleStep = 1.0;  // disable safety clamp for unit test
  implicitParams.fluidTorqueScale = 0.0;
  std::vector<T> zeroFluidVec(config.nSegments, 0.0);
  for (int s = 0; s < 20; ++s) {
    advanceSoftBackboneImplicit(config, restingPreferredState, zeroFluidVec,
                                1e-4, implicitParams, implicitState);
  }
  T maxRestDrift = 0.0;
  for (int i = 0; i < config.nSegments; ++i) {
    maxRestDrift = std::max(maxRestDrift, std::abs(implicitState.theta[i]));
    maxRestDrift = std::max(maxRestDrift, std::abs(implicitState.omega[i]));
  }
  assert(maxRestDrift < 1e-12);

  // (b) With non-trivial preferred and zero fluid, segment 0 must be pinned
  //     exactly to preferred at every step.
  SoftBackboneState implicitFromStraight = straight;
  advanceSoftBackboneImplicit(config, preferredState, zeroFluidVec,
                              1e-4, implicitParams, implicitFromStraight);
  assert(std::abs(implicitFromStraight.theta.front() -
                  preferredState.theta.front()) < 1e-12);
  assert(std::abs(implicitFromStraight.omega.front() -
                  preferredState.omega.front()) < 1e-12);

  // (c) Inertia / phase-lag check: relax a perturbation at zero preferred
  //     with no fluid load — energy must not grow over a short integration
  //     (numerical dissipation of implicit Euler is harmless here).
  SoftBackboneState perturbed = restingPreferredState;
  perturbed.theta[config.nSegments / 2] = 0.05;  // small bend
  T initialMaxAbs = 0.05;
  for (int s = 0; s < 50; ++s) {
    advanceSoftBackboneImplicit(config, restingPreferredState, zeroFluidVec,
                                1e-4, implicitParams, perturbed);
  }
  T finalMaxAbs = 0.0;
  for (T th : perturbed.theta) {
    finalMaxAbs = std::max(finalMaxAbs, std::abs(th));
  }
  assert(finalMaxAbs <= initialMaxAbs + 1e-9);

  // (d) Fluid torque actually drives joint angle changes when scale > 0.
  SoftBackboneState forced = preferredState;
  SoftBackboneDynamicsParams forcedParams = implicitParams;
  forcedParams.fluidTorqueScale = 1.0;
  std::vector<T> fluidLoad(config.nSegments, 0.0);
  fluidLoad[config.nSegments - 2] = 0.001;
  fluidLoad[config.nSegments - 1] = -0.001;
  const SoftBackboneDynamicsDiagnostics forcedDiag =
    advanceSoftBackboneImplicit(config, preferredState, fluidLoad,
                                1e-4, forcedParams, forced);
  assert(forcedDiag.maxAbsFluidSegmentTorqueNm > 0.0);
  assert(forcedDiag.maxAbsAngleStep > 0.0);
}
