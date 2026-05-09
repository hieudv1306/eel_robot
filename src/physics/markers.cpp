#include "physics/markers.hpp"

#include "physics/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

int LagrangianMarkers::size() const { return static_cast<int>(x.size()); }

void LagrangianMarkers::resize(int n) {
  x.resize(n); y.resize(n);
  ud.resize(n); vd.resize(n);
  ds.resize(n);
  fx.resize(n); fy.resize(n);
}

void buildLagrangianMarkers(
    const EelParams& p,
    T t, T Vx, T Vy, T omegaZ,
    T xCm, T yCm, T theta,
    T dtLbm, T ampRamp,
    LagrangianMarkers& markers)
{
  std::vector<T> globX, globY, vxDef, vyDef, allDs;
  buildCapsuleGeometry(p, t, xCm, yCm, theta, dtLbm, ampRamp, globX, globY, vxDef, vyDef, allDs);

  int n = static_cast<int>(globX.size());
  markers.resize(n);

  for (int i = 0; i < n; ++i) {
    markers.x[i]  = globX[i];
    markers.y[i]  = globY[i];
    markers.ds[i] = allDs[i];

    // Desired velocity = rigid-body + deformation.
    // Rigid-body state is advanced in lattice-step units (no extra dt factor).
    T rx = globX[i] - xCm;
    T ry = globY[i] - yCm;
    markers.ud[i] = Vx - omegaZ * ry + vxDef[i];
    markers.vd[i] = Vy + omegaZ * rx + vyDef[i];
  }
}

void buildLagrangianMarkersFromSoftBackbone(
    const EelParams& p,
    const SoftBackboneConfig& config,
    T t, T Vx, T Vy, T omegaZ,
    T xCm, T yCm, T theta,
    T dtLbm, T ampRamp,
    LagrangianMarkers& markers)
{
  const SoftBackboneState state =
    preferredBackboneStateWave(p, config, t, dtLbm, ampRamp);
  const SoftBackboneState statePlus =
    preferredBackboneStateWave(p, config, t + dtLbm, dtLbm, ampRamp);
  const SoftBackboneState stateMinus =
    preferredBackboneStateWave(p, config, t - dtLbm, dtLbm, ampRamp);

  std::vector<T> globX, globY, vxDef, vyDef, allDs;
  buildCapsuleGeometryFromBackboneMotion(
    p, config, state, statePlus, stateMinus,
    xCm, yCm, theta, globX, globY, vxDef, vyDef, allDs);

  const int n = static_cast<int>(globX.size());
  markers.resize(n);

  for (int i = 0; i < n; ++i) {
    markers.x[i]  = globX[i];
    markers.y[i]  = globY[i];
    markers.ds[i] = allDs[i];

    T rx = globX[i] - xCm;
    T ry = globY[i] - yCm;
    markers.ud[i] = Vx - omegaZ * ry + vxDef[i];
    markers.vd[i] = Vy + omegaZ * rx + vyDef[i];
  }
}

void buildLagrangianMarkersFromSoftBackboneState(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    T Vx, T Vy, T omegaZ,
    T xCm, T yCm, T theta,
    T dtLbm,
    LagrangianMarkers& markers)
{
  const SoftBackboneState statePlus =
    extrapolateBackboneState(state, dtLbm);
  const SoftBackboneState stateMinus =
    extrapolateBackboneState(state, -dtLbm);

  std::vector<T> globX, globY, vxDef, vyDef, allDs;
  buildCapsuleGeometryFromBackboneMotion(
    p, config, state, statePlus, stateMinus,
    xCm, yCm, theta, globX, globY, vxDef, vyDef, allDs);

  const int n = static_cast<int>(globX.size());
  markers.resize(n);

  for (int i = 0; i < n; ++i) {
    markers.x[i]  = globX[i];
    markers.y[i]  = globY[i];
    markers.ds[i] = allDs[i];

    T rx = globX[i] - xCm;
    T ry = globY[i] - yCm;
    markers.ud[i] = Vx - omegaZ * ry + vxDef[i];
    markers.vd[i] = Vy + omegaZ * rx + vyDef[i];
  }
}

SoftBackboneForceProjection projectMarkerForcesToSoftBackbone(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    T xCm, T yCm, T theta,
    const LagrangianMarkers& markers)
{
  SoftBackboneForceProjection out;
  if (!config.valid || config.nSegments <= 0 || markers.size() <= 0) {
    return out;
  }

  const SoftBackboneCenterline centerline =
    buildSoftBackboneCenterlineLU(p, config, state, xCm, yCm, theta);
  const int nSegments = static_cast<int>(centerline.segmentX.size());
  if (nSegments <= 0 ||
      static_cast<int>(centerline.nodeX.size()) != nSegments + 1) {
    return out;
  }

  out.segmentTorqueNm.assign(nSegments, T(0));
  std::vector<T> segmentTorqueLat(nSegments, T(0));

  const T dxM = (p.totalGeometricLengthLU() > T(0))
              ? p.physicalBodyLengthM / p.totalGeometricLengthLU() : T(0);
  const T dtSec = p.dtLbm();
  const T thicknessM = (p.bodyThicknessM > T(0))
                     ? p.bodyThicknessM : dxM;
  const T forceScaleN =
    (dxM > T(0) && dtSec > T(0) && p.fluidDensityKgM3 > T(0))
      ? p.fluidDensityKgM3 * thicknessM * dxM * dxM * dxM
        / (dtSec * dtSec)
      : T(0);
  out.latticeTorqueToNm = forceScaleN * dxM;

  for (int m = 0; m < markers.size(); ++m) {
    const T xm = markers.x[m];
    const T ym = markers.y[m];
    int bestSegment = 0;
    T bestDist2 = std::numeric_limits<T>::max();
    for (int s = 0; s < nSegments; ++s) {
      const T ax = centerline.nodeX[s];
      const T ay = centerline.nodeY[s];
      const T bx = centerline.nodeX[s + 1];
      const T by = centerline.nodeY[s + 1];
      const T vx = bx - ax;
      const T vy = by - ay;
      const T denom = vx * vx + vy * vy;
      T u = T(0);
      if (denom > T(1e-12)) {
        u = ((xm - ax) * vx + (ym - ay) * vy) / denom;
        u = std::clamp(u, T(0), T(1));
      }
      const T px = ax + u * vx;
      const T py = ay + u * vy;
      const T dx = xm - px;
      const T dy = ym - py;
      const T dist2 = dx * dx + dy * dy;
      if (dist2 < bestDist2) {
        bestDist2 = dist2;
        bestSegment = s;
      }
    }

    // markers.fx/fy are forces applied to the fluid. The body reaction is the
    // opposite sign, integrated with the marker arc-length weight.
    const T bodyFxLat = -markers.fx[m] * markers.ds[m];
    const T bodyFyLat = -markers.fy[m] * markers.ds[m];
    out.netForceXLat += bodyFxLat;
    out.netForceYLat += bodyFyLat;
    const T rx = xm - centerline.segmentX[bestSegment];
    const T ry = ym - centerline.segmentY[bestSegment];
    const T torqueLat = rx * bodyFyLat - ry * bodyFxLat;
    segmentTorqueLat[bestSegment] += torqueLat;
    out.netTorqueLat += torqueLat;
  }

  for (int s = 0; s < nSegments; ++s) {
    out.segmentTorqueNm[s] = segmentTorqueLat[s] * out.latticeTorqueToNm;
    out.maxAbsSegmentTorqueNm =
      std::max(out.maxAbsSegmentTorqueNm, std::abs(out.segmentTorqueNm[s]));
  }
  return out;
}
