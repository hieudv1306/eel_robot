#include "physics/markers.hpp"

#include "physics/geometry.hpp"

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
