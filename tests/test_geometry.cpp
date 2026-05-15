#include "core/params.hpp"
#include "physics/geometry.hpp"
#include "physics/markers.hpp"
#include "physics/material.hpp"
#include "physics/soft_rod.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

int main() {
  EelParams p;
  p.nSpine = 20;
  p.bodyRadius = 3.0;
  p.eelScale = 40.0;
  p.useAspectRatioGeometry = false;

  std::vector<T> x, y, ux, uy, ds;
  buildCapsuleGeometry(p, 0.25, 100.0, 60.0, 0.1, p.dtLbm(), 1.0,
                       x, y, ux, uy, ds);

  assert(!x.empty());
  assert(x.size() == y.size());
  assert(x.size() == ux.size());
  assert(x.size() == uy.size());
  assert(x.size() == ds.size());

  std::vector<T> xp, yp, uxp, uyp, dsp;
  std::vector<T> xm, ym, uxm, uym, dsm;
  const T t = 1.0;
  const T dt = p.dtLbm();
  buildCapsuleGeometry(p, t, 100.0, 60.0, 0.0, dt, 1.0,
                       x, y, ux, uy, ds);
  buildCapsuleGeometry(p, t + dt, 100.0, 60.0, 0.0, dt, 1.0,
                       xp, yp, uxp, uyp, dsp);
  buildCapsuleGeometry(p, t - dt, 100.0, 60.0, 0.0, dt, 1.0,
                       xm, ym, uxm, uym, dsm);
  assert(x.size() == xp.size());
  assert(x.size() == xm.size());
  T maxAbsUx = 0.0;
  T maxAbsUy = 0.0;
  for (size_t i = 0; i < x.size(); ++i) {
    const T uxFd = 0.5 * (xp[i] - xm[i]);
    const T uyFd = 0.5 * (yp[i] - ym[i]);
    assert(std::abs(ux[i] - uxFd) < 1e-12);
    assert(std::abs(uy[i] - uyFd) < 1e-12);
    maxAbsUx = std::max(maxAbsUx, std::abs(ux[i]));
    maxAbsUy = std::max(maxAbsUy, std::abs(uy[i]));
  }
  assert(maxAbsUx > 1e-8);
  assert(maxAbsUy > 1e-8);

  EelParams pInext = p;
  pInext.nSpine = 50;
  buildCapsuleGeometry(pInext, 2.0, 100.0, 60.0, 0.0, dt, 1.0,
                       x, y, ux, uy, ds);
  const T dsSpine = pInext.eelScale * pInext.eelLength / (pInext.nSpine - 1);
  const int nTailCap = std::max(
    static_cast<int>(std::ceil(3.14159265358979323846 * pInext.bodyRadius / dsSpine)),
    8);
  const int lowerStart = pInext.nSpine + (nTailCap - 1);
  std::vector<T> cx(pInext.nSpine), cy(pInext.nSpine);
  for (int i = 0; i < pInext.nSpine; ++i) {
    const int lowerIdx = lowerStart + (pInext.nSpine - 1 - i);
    cx[i] = T(0.5) * (x[i] + x[lowerIdx]);
    cy[i] = T(0.5) * (y[i] + y[lowerIdx]);
  }
  T centerlineArc = 0.0;
  for (int i = 1; i < pInext.nSpine; ++i) {
    const T dx = cx[i] - cx[i - 1];
    const T dy = cy[i] - cy[i - 1];
    centerlineArc += std::sqrt(dx * dx + dy * dy);
  }
  assert(std::abs(centerlineArc - pInext.centerlineLengthLU()) < 1e-10);

  auto e = bodyForwardAxis(0.3);
  auto n = bodyLateralAxis(0.3);
  const T dot = e[0] * n[0] + e[1] * n[1];
  const T eNorm = std::sqrt(e[0] * e[0] + e[1] * e[1]);
  const T nNorm = std::sqrt(n[0] * n[0] + n[1] * n[1]);
  assert(std::abs(dot) < 1e-12);
  assert(std::abs(eNorm - 1.0) < 1e-12);
  assert(std::abs(nNorm - 1.0) < 1e-12);

  pInext.physicalBodyLengthM = 0.30;
  pInext.bodyThicknessM = 0.02;
  const MaterialProperties material = resolveMaterialProperties(pInext);
  const PlanarRodSectionEstimate rod =
    estimatePlanarRodSection(pInext, material);
  const SoftBackboneConfig backbone =
    makeSoftBackboneConfig(pInext, rod, pInext.nSpine - 1);
  assert(backbone.valid);
  const SoftBackboneState preferred =
    preferredBackboneStateWave(pInext, backbone, 2.0, dt);

  EelParams pReverse = pInext;
  pReverse.waveDirection = WaveDirection::TailToHead;
  const SoftBackboneState reversed =
    preferredBackboneStateWave(pReverse, backbone, 2.0, dt);
  assert(reversed.theta.size() == preferred.theta.size());
  T maxWaveDirectionThetaDiff = 0.0;
  for (size_t i = 0; i < preferred.theta.size(); ++i) {
    maxWaveDirectionThetaDiff =
      std::max(maxWaveDirectionThetaDiff,
               std::abs(preferred.theta[i] - reversed.theta[i]));
  }
  assert(maxWaveDirectionThetaDiff > 1e-8);

  std::vector<T> xb, yb, dsb;
  buildCapsuleGeometryFromBackboneState(
    pInext, backbone, preferred, 100.0, 60.0, 0.0, xb, yb, dsb);
  assert(!xb.empty());
  assert(xb.size() == yb.size());
  assert(xb.size() == dsb.size());

  std::vector<T> xpB, ypB, uxB, uyB, dsMotionB;
  const SoftBackboneState preferredPlus =
    preferredBackboneStateWave(pInext, backbone, 2.0 + dt, dt);
  const SoftBackboneState preferredMinus =
    preferredBackboneStateWave(pInext, backbone, 2.0 - dt, dt);
  buildCapsuleGeometryFromBackboneMotion(
    pInext, backbone, preferred, preferredPlus, preferredMinus,
    100.0, 60.0, 0.0, xpB, ypB, uxB, uyB, dsMotionB);
  assert(xpB.size() == xb.size());
  assert(xpB.size() == uxB.size());
  T maxSoftSpeed = 0.0;
  for (size_t i = 0; i < uxB.size(); ++i) {
    maxSoftSpeed = std::max(maxSoftSpeed, std::sqrt(uxB[i] * uxB[i] +
                                                    uyB[i] * uyB[i]));
  }
  assert(maxSoftSpeed > 1e-8);

  LagrangianMarkers markers;
  buildLagrangianMarkersFromSoftBackboneState(
    pInext, backbone, preferred, 0.0, 0.0, 0.0,
    100.0, 60.0, 0.0, dt, markers);
  assert(markers.size() > 0);
  for (int i = 0; i < markers.size(); ++i) {
    markers.fx[i] = 0.0;
    markers.fy[i] = 0.0;
  }
  auto noLoad = projectMarkerForcesToSoftBackbone(
    pInext, backbone, preferred, 100.0, 60.0, 0.0, markers);
  assert(noLoad.segmentTorqueNm.size() == static_cast<size_t>(backbone.nSegments));
  assert(noLoad.maxAbsSegmentTorqueNm == 0.0);

  markers.fx[0] = 1.0;
  markers.fy[0] = -0.5;
  auto loaded = projectMarkerForcesToSoftBackbone(
    pInext, backbone, preferred, 100.0, 60.0, 0.0, markers);
  assert(loaded.segmentTorqueNm.size() == static_cast<size_t>(backbone.nSegments));
  assert(loaded.latticeTorqueToNm > 0.0);
  assert(loaded.maxAbsSegmentTorqueNm > 0.0);

  auto virtualLoad = projectMarkerForcesToSoftBackbone(
    pInext, backbone, preferred, 100.0, 60.0, 0.0, markers,
    SoftBackboneLoadProjection::CrossSectionVirtualWork);
  assert(virtualLoad.segmentTorqueNm.size() ==
         static_cast<size_t>(backbone.nSegments));
  assert(virtualLoad.latticeTorqueToNm > 0.0);
  assert(virtualLoad.maxAbsSegmentTorqueNm > 0.0);
  assert(std::abs(virtualLoad.maxAbsSegmentTorqueNm -
                  loaded.maxAbsSegmentTorqueNm) > 1e-18);

  // ----- Soft transitional path matches legacy prescribed-wave -----
  // When --bodyKinematics=soft_backbone --softBackboneDynamics=false the
  // backbone state is the analytical preferred wave, so the IBM markers
  // produced by the soft path must agree with the legacy prescribed-wave
  // markers to O(ds^2).  A larger gap means the segment-angle -> centerline
  // conversion has lost boundary tangent accuracy (issue 1 in the
  // 2026-05-09 validation).
  EelParams pCmp = pInext;
  pCmp.nSpine = 100;
  pCmp.bodyRadius = 4.0;
  pCmp.eelScale = 60.0;
  pCmp.eelFreq = 1.6;
  pCmp.eelLambda = 1.0;
  pCmp.eelA0 = 0.07;
  pCmp.restTime = 0.1;
  pCmp.rampTime = 0.5;
  pCmp.dtAnim = 0.04;
  pCmp.substeps = 20;
  const T tCmp = 1.0;
  const T dtCmp = pCmp.dtLbm();
  const MaterialProperties matCmp = resolveMaterialProperties(pCmp);
  const PlanarRodSectionEstimate rodCmp = estimatePlanarRodSection(pCmp, matCmp);
  const SoftBackboneConfig backboneCmp =
    makeSoftBackboneConfig(pCmp, rodCmp, pCmp.nSpine - 1);
  assert(backboneCmp.valid);

  LagrangianMarkers legacyMarkers, softMarkers;
  buildLagrangianMarkers(pCmp, tCmp, 0, 0, 0, 100.0, 60.0, 0.0,
                         dtCmp, 1.0, legacyMarkers);
  buildLagrangianMarkersFromSoftBackbone(
    pCmp, backboneCmp, tCmp, 0, 0, 0, 100.0, 60.0, 0.0,
    dtCmp, 1.0, softMarkers);
  assert(legacyMarkers.size() == softMarkers.size());

  T maxDxAbs = 0.0, maxDyAbs = 0.0;
  T maxDuAbs = 0.0, maxDvAbs = 0.0;
  for (int i = 0; i < legacyMarkers.size(); ++i) {
    maxDxAbs = std::max(maxDxAbs,
                        std::abs(softMarkers.x[i] - legacyMarkers.x[i]));
    maxDyAbs = std::max(maxDyAbs,
                        std::abs(softMarkers.y[i] - legacyMarkers.y[i]));
    maxDuAbs = std::max(maxDuAbs,
                        std::abs(softMarkers.ud[i] - legacyMarkers.ud[i]));
    maxDvAbs = std::max(maxDvAbs,
                        std::abs(softMarkers.vd[i] - legacyMarkers.vd[i]));
  }
  // Tolerances chosen with margin around the residual O(ds^2) discretisation
  // gap (~5e-3 LU positions, ~3e-4 LU/step velocities at nSpine=100).  If a
  // future change reverts the segment->node tangent conversion, these will
  // immediately blow past 1e-1 LU and catch the regression.
  assert(maxDxAbs < 1e-2);
  assert(maxDyAbs < 1e-2);
  assert(maxDuAbs < 1e-3);
  assert(maxDvAbs < 1e-3);
}
