#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>
#include <vector>

#include "core/params.hpp"
#include "core/types.hpp"
#include "physics/markers.hpp"
#include "solver/openlb_setup.hpp"

// ============================================================
//  Peskin 4-point delta function
// ============================================================
inline T peskinDelta(T r) {
  T ar = std::abs(r);
  if (ar < 1.0) {
    return 0.125 * (3.0 - 2.0 * ar + std::sqrt(1.0 + 4.0 * ar - 4.0 * ar * ar));
  } else if (ar < 2.0) {
    return 0.125 * (5.0 - 2.0 * ar - std::sqrt(-7.0 + 12.0 * ar - 4.0 * ar * ar));
  }
  return 0.0;
}

inline T delta2D(T dx, T dy) {
  return peskinDelta(dx) * peskinDelta(dy);
}

// ============================================================
//  IBM operations on OpenLB lattice fields (lattice-units)
// ============================================================
struct IbmResult {
  T forceX  = 0.0;
  T forceY  = 0.0;
  T torque  = 0.0;
  // power is the *constraint* power F·U_marker (legacy semantics).  Renamed
  // logically to constraintPower in v6, but the field is left as `power` for
  // backward compatibility with all the downstream accumulators that already
  // read it.  The decomposition powerRigid + powerDef ≈ power holds
  // exactly in lattice arithmetic because U_marker = U_rigid + U_def.
  T power       = 0.0;  // F · U_marker  (constraint power, legacy)
  T powerRigid  = 0.0;  // F · U_rigid   (rigid-body component)
  T powerDef    = 0.0;  // F · (U_marker - U_rigid) = F · U_def
  // Slip statistics are unweighted per-marker averages of the IBM velocity
  // mismatch |U_desired-U_interp|. They quantify enforcement quality only.
  T meanSlipMag = 0.0;
  T maxSlipMag = 0.0;
  // Residual slip = mismatch left over after the multi-direct-forcing
  // iterations, measured against a marker velocity re-interpolated from the
  // sparse Eulerian correction buffer.  This is an IBM-iteration convergence
  // metric; it is NOT the post-fluid-update slip (which only becomes
  // observable on the next IBM pass).
  T meanResidualSlipMag = 0.0;
  T maxResidualSlipMag  = 0.0;
  // Marker-force statistics are unweighted per-marker averages of the direct-
  // forcing IBM line-force density magnitude |f_marker| before ds integration.
  // They are not total body-force proxies; the integrated reaction force uses
  // bodyFx/bodyFy below after multiplication by the marker arc-length weight ds.
  T meanMarkerForceMag = 0.0;
  T maxMarkerForceMag = 0.0;
  T meanDesiredMarkerSpeedMag = 0.0;
  T normalizedMeanSlip = 0.0;
};

// Interpolate velocity at marker position from the OpenLB lattice, compute
// the multi-direct-forcing (MDF) IBM correction, and spread the accumulated
// force into the OpenLB FORCE field.
//
//  Forcing law (v6, physicalized):
//    f_marker = alphaIBM * rho_local * (U_desired - U_estimate) / dt_effective
//  Numerical assumptions used here:
//    * rho_local = 1.0 — the OpenLB lattice runs near-incompressible with
//      ρ ≈ 1 in lattice units (the unit-convention block at the top of the
//      file documents this).  No per-marker density interpolation is done so
//      the interpolation/spreading stencils remain unchanged.
//    * dt_effective = 1.0 — Newton-Euler and gait integration in this file
//      operate with δt = 1 lattice step throughout.  The "physicalized"
//      coefficient alphaIBM therefore equals the classical penalty gain in
//      the non-dimensional setting (alphaIBM=2 ≡ legacy --kappa=2).
//
//  Multi direct forcing iterations (Wang & Zhang 2011 style, Eulerian-buffer
//  approximation):
//    * One velocity interpolation per marker per IBM call.
//    * Per-iteration: residual = U_desired - U_estimate, dF accumulates.
//      The incremental force is spread to a sparse Eulerian velocity-
//      correction buffer and subsequent iterations re-interpolate from that
//      buffer.  This captures marker cross-talk through the same Peskin kernel
//      used for the final OpenLB FORCE field, without advancing the fluid.
//    * The accumulated marker force is spread into the FORCE field once,
//      after the loop.  This preserves the existing
//      single-IBM-per-collideAndStream coupling order.
template <typename ResetForceField>
void ibmStep(
    const EelParams& p,
    LagrangianMarkers& markers,
    T xCm, T yCm,
    T Vx, T Vy, T omegaZ,
    olb::SuperLattice<T, DESCRIPTOR>& sLattice,
    T alphaIBM,
    int ibmIterations,
    ResetForceField&& resetForceField,
    IbmResult& result)
{
  const int nx = p.nx;
  const int ny = p.ny;
  const int nM = markers.size();
  const int nIters = std::max(1, ibmIterations);

  result = {};
  sLattice.setProcessingContext(olb::ProcessingContext::Evaluation);

  // Numerical assumptions for the physicalized IBM (see function header).
  const T rhoLocal = T(1);
  const T dtEff    = T(1);
  const T gain     = alphaIBM * rhoLocal / dtEff;
  const T kick     = dtEff / rhoLocal;  // ΔU per unit Eulerian force density

  struct MarkerStencil {
    int n = 0;
    std::array<int, 16> cell{};
    std::array<T, 16> weight{};
  };
  struct VelocityCorrection {
    T u = T(0);
    T v = T(0);
  };

  // First pass: cache the Peskin support and interpolate the fluid velocity at
  // every marker exactly once.  Later IBM iterations re-interpolate only the
  // sparse Eulerian correction field.
  std::vector<T> uInterp(nM, T(0));
  std::vector<T> vInterp(nM, T(0));
  std::vector<T> uEst(nM, T(0));
  std::vector<T> vEst(nM, T(0));
  std::vector<T> fxAcc(nM, T(0));
  std::vector<T> fyAcc(nM, T(0));
  std::vector<T> dFx(nM, T(0));
  std::vector<T> dFy(nM, T(0));
  std::vector<int> i0Marker(nM, 0);
  std::vector<int> j0Marker(nM, 0);
  std::vector<MarkerStencil> stencils(nM);

  for (int m = 0; m < nM; ++m) {
    const T xm = markers.x[m];
    const T ym = markers.y[m];
    const int i0 = static_cast<int>(std::floor(xm)) - 1;
    const int j0 = static_cast<int>(std::floor(ym)) - 1;
    i0Marker[m] = i0;
    j0Marker[m] = j0;

    T uI = T(0), vI = T(0);
    for (int di = 0; di < 4; ++di) {
      for (int dj = 0; dj < 4; ++dj) {
        const int ii = i0 + di;
        const int jj = j0 + dj;
        if (ii >= 0 && ii < nx && jj >= 0 && jj < ny) {
          const T d = delta2D(xm - T(ii), ym - T(jj));
          MarkerStencil& stencil = stencils[m];
          stencil.cell[stencil.n] = ii * ny + jj;
          stencil.weight[stencil.n] = d;
          ++stencil.n;
          T uCell[2] = {0.0, 0.0};
          sLattice.get(0, ii, jj).computeU(uCell);
          uI += uCell[0] * d;
          vI += uCell[1] * d;
        }
      }
    }
    uInterp[m] = uI;
    vInterp[m] = vI;
    uEst[m] = uI;
    vEst[m] = vI;
  }

  std::unordered_map<int, VelocityCorrection> eulerianCorrection;
  eulerianCorrection.reserve(static_cast<size_t>(std::max(1, nM) * 8));

  auto interpolateCorrection = [&](int m, T& uCorr, T& vCorr) {
    uCorr = T(0);
    vCorr = T(0);
    const MarkerStencil& stencil = stencils[m];
    for (int k = 0; k < stencil.n; ++k) {
      const auto it = eulerianCorrection.find(stencil.cell[k]);
      if (it == eulerianCorrection.end()) continue;
      const T w = stencil.weight[k];
      uCorr += it->second.u * w;
      vCorr += it->second.v * w;
    }
  };

  auto spreadCorrection = [&](int m, T fx, T fy) {
    const MarkerStencil& stencil = stencils[m];
    const T dsM = markers.ds[m];
    for (int k = 0; k < stencil.n; ++k) {
      const T scaledWeight = dsM * stencil.weight[k] * kick;
      auto& corr = eulerianCorrection[stencil.cell[k]];
      corr.u += fx * scaledWeight;
      corr.v += fy * scaledWeight;
    }
  };

  // MDF iterations.  N=1 still reproduces the classical single-pass direct
  // forcing exactly in the final FORCE field.  N>1 now re-spreads correction
  // increments to a sparse Eulerian buffer and re-interpolates from it, so
  // marker cross-coupling through the Peskin kernel is represented.
  for (int iter = 0; iter < nIters; ++iter) {
    for (int m = 0; m < nM; ++m) {
      T uCorr = T(0), vCorr = T(0);
      interpolateCorrection(m, uCorr, vCorr);
      uEst[m] = uInterp[m] + uCorr;
      vEst[m] = vInterp[m] + vCorr;
      const T sx = markers.ud[m] - uEst[m];
      const T sy = markers.vd[m] - vEst[m];
      dFx[m] = gain * sx;
      dFy[m] = gain * sy;
    }
    for (int m = 0; m < nM; ++m) {
      fxAcc[m] += dFx[m];
      fyAcc[m] += dFy[m];
      spreadCorrection(m, dFx[m], dFy[m]);
    }
  }

  for (int m = 0; m < nM; ++m) {
    T uCorr = T(0), vCorr = T(0);
    interpolateCorrection(m, uCorr, vCorr);
    uEst[m] = uInterp[m] + uCorr;
    vEst[m] = vInterp[m] + vCorr;
  }

  // Aggregate per-marker statistics, body-frame reaction force/torque, and
  // the power decomposition.  Power is integrated over markers with the
  // arc-length weight ds, matching the legacy summation.
  T slipSum = T(0);
  T residualSlipSum = T(0);
  T markerForceSum  = T(0);
  T desiredSpeedSum = T(0);

  for (int m = 0; m < nM; ++m) {
    const T fxM = fxAcc[m];
    const T fyM = fyAcc[m];
    markers.fx[m] = fxM;
    markers.fy[m] = fyM;

    const T slipX0 = markers.ud[m] - uInterp[m];
    const T slipY0 = markers.vd[m] - vInterp[m];
    const T slipMag = std::sqrt(slipX0 * slipX0 + slipY0 * slipY0);

    const T resX = markers.ud[m] - uEst[m];
    const T resY = markers.vd[m] - vEst[m];
    const T residualMag = std::sqrt(resX * resX + resY * resY);

    const T markerForceMag = std::sqrt(fxM * fxM + fyM * fyM);
    const T desiredSpeedMag = std::sqrt(markers.ud[m] * markers.ud[m]
                                      + markers.vd[m] * markers.vd[m]);

    slipSum         += slipMag;
    residualSlipSum += residualMag;
    markerForceSum  += markerForceMag;
    desiredSpeedSum += desiredSpeedMag;
    result.maxSlipMag         = std::max(result.maxSlipMag, slipMag);
    result.maxResidualSlipMag = std::max(result.maxResidualSlipMag, residualMag);
    result.maxMarkerForceMag  = std::max(result.maxMarkerForceMag, markerForceMag);

    const T xm = markers.x[m];
    const T ym = markers.y[m];
    const T dsM = markers.ds[m];

    // Body reaction (opposite sign, integrated with arc-length weight).
    const T bodyFx = -fxM * dsM;
    const T bodyFy = -fyM * dsM;
    result.forceX += bodyFx;
    result.forceY += bodyFy;

    const T rx = xm - xCm;
    const T ry = ym - yCm;
    result.torque += rx * bodyFy - ry * bodyFx;

    // Power decomposition: U_marker = U_rigid + U_def.
    //   U_rigid = (Vx - omegaZ·ry, Vy + omegaZ·rx)   (centre-of-mass
    //              translation + rigid rotation about CM)
    //   U_def   = U_marker - U_rigid                 (deformation only)
    // The identity P_total = P_rigid + P_def holds exactly in lattice
    // arithmetic because the decomposition is purely additive.
    const T uRigidX = Vx - omegaZ * ry;
    const T uRigidY = Vy + omegaZ * rx;
    const T uDefX   = markers.ud[m] - uRigidX;
    const T uDefY   = markers.vd[m] - uRigidY;

    result.power      += (fxM * markers.ud[m] + fyM * markers.vd[m]) * dsM;
    result.powerRigid += (fxM * uRigidX       + fyM * uRigidY      ) * dsM;
    result.powerDef   += (fxM * uDefX         + fyM * uDefY        ) * dsM;
  }

  if (nM > 0) {
    const T invMarkers = T(1) / T(nM);
    result.meanSlipMag         = slipSum         * invMarkers;
    result.meanResidualSlipMag = residualSlipSum * invMarkers;
    result.meanMarkerForceMag  = markerForceSum  * invMarkers;
    result.meanDesiredMarkerSpeedMag = desiredSpeedSum * invMarkers;
    if (result.meanDesiredMarkerSpeedMag > T(1e-12)) {
      result.normalizedMeanSlip =
        result.meanSlipMag / result.meanDesiredMarkerSpeedMag;
    }
  }

  // Replace the previous IBM Eulerian force field with the newly computed
  // marker force.  The time loop intentionally performs one IBM evaluation
  // per collideAndStream(); the MDF iterations above happen in a sparse
  // Eulerian correction buffer, and the final accumulated marker force is
  // spread into OpenLB once.
  resetForceField();

  for (int m = 0; m < nM; ++m) {
    const T xm = markers.x[m];
    const T ym = markers.y[m];
    const int i0 = i0Marker[m];
    const int j0 = j0Marker[m];
    const T fxM = markers.fx[m];
    const T fyM = markers.fy[m];
    const T dsM = markers.ds[m];

    // Spread force into the OpenLB lattice FORCE field (Peskin 4-point
    // delta — kernel unchanged from v5).
    for (int di = 0; di < 4; ++di) {
      for (int dj = 0; dj < 4; ++dj) {
        const int ii = i0 + di;
        const int jj = j0 + dj;
        if (ii >= 0 && ii < nx && jj >= 0 && jj < ny) {
          const T d = delta2D(xm - T(ii), ym - T(jj));
          auto cell = sLattice.get(0, ii, jj);
          auto force = cell.template getField<olb::descriptors::FORCE>();
          force[0] += fxM * dsM * d;
          force[1] += fyM * dsM * d;
          cell.template setField<olb::descriptors::FORCE>(force);
        }
      }
    }
  }
  sLattice.template setProcessingContext<
    olb::Array<olb::descriptors::FORCE>>(olb::ProcessingContext::Simulation);
}
