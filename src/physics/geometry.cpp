#include "physics/geometry.hpp"

#include "physics/gait.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

// ============================================================
//  Build capsule geometry and lagrangian markers
// ============================================================
namespace {

void buildCapsulePositionFrame(
    const EelParams& p,
    T t, T xCm, T yCm, T theta,
    T ampRamp,
    std::vector<T>& globX, std::vector<T>& globY,
    std::vector<T>& allDs)
{
  const T L      = p.eelLength;
  const int nSp  = p.nSpine;
  const T scale  = p.eelScale;
  const T freq   = p.eelFreq;
  const T lam    = p.eelLambda;
  const T A0     = p.eelA0;
  const T r      = p.bodyRadius;
  const T pi2    = 2.0 * M_PI;

  // Spine parameter: s=0 at the head, s=L at the tail.
  std::vector<T> s(nSp);
  for (int i = 0; i < nSp; ++i) {
    s[i] = L * i / (nSp - 1);
  }

  auto [tActive, gaitRamp, gaitRampRate] =
    actuationProfile(t, p.restTime, p.rampTime);
  (void)gaitRampRate;
  const T ampScale = ampRamp * gaitRamp;

  // Local frame: centered at CM (s = L/2)
  const T sCm = 0.5 * L;

  std::vector<T> cxL(nSp), cyL(nSp);
  std::vector<T> Abase(nSp), Aval(nSp), phase(nSp);
  std::vector<T> txV(nSp), tyV(nSp), nxN(nSp), nyN(nSp);

  for (int i = 0; i < nSp; ++i) {
    // s/lambda - f t sends the body wave from head to tail (increasing s).
    phase[i] = pi2 * (s[i] / lam - freq * tActive);
    Abase[i] = amplitudeEnvelope(s[i], L, A0);
    Aval[i]  = Abase[i] * ampScale;
  }

  const T dsParam = L / (nSp - 1);
  const T dsSpine = scale * L / (nSp - 1);

  auto amplitudeDerivative = [&](int i) -> T {
    if (i > 0 && i < nSp - 1) {
      return (Aval[i+1] - Aval[i-1]) / (T(2) * dsParam);
    }
    if (i == 0) {
      return (Aval[1] - Aval[0]) / dsParam;
    }
    return (Aval[nSp-1] - Aval[nSp-2]) / dsParam;
  };

  std::vector<T> slope(nSp);
  for (int i = 0; i < nSp; ++i) {
    const T dAds = amplitudeDerivative(i);
    slope[i] =
      Aval[i] * std::cos(phase[i]) * (pi2 / lam)
      + dAds * std::sin(phase[i]);
  }

  if (p.geometryKinematics == GeometryKinematics::InextensibleWave) {
    std::vector<T> tangentAngle(nSp);
    for (int i = 0; i < nSp; ++i) {
      tangentAngle[i] = std::atan(slope[i]);
      txV[i] = std::cos(tangentAngle[i]);
      tyV[i] = std::sin(tangentAngle[i]);
      nxN[i] = -tyV[i];
      nyN[i] =  txV[i];
    }

    cxL[0] = T(0);
    cyL[0] = T(0);
    for (int i = 1; i < nSp; ++i) {
      const T aMid = T(0.5) * (tangentAngle[i - 1] + tangentAngle[i]);
      cxL[i] = cxL[i - 1] + dsSpine * std::cos(aMid);
      cyL[i] = cyL[i - 1] + dsSpine * std::sin(aMid);
    }

    const T midIndex = T(0.5) * T(nSp - 1);
    const int iMid0 = static_cast<int>(std::floor(midIndex));
    const int iMid1 = std::min(iMid0 + 1, nSp - 1);
    const T midFrac = midIndex - T(iMid0);
    const T xMid = cxL[iMid0] * (T(1) - midFrac) + cxL[iMid1] * midFrac;
    const T yMid = cyL[iMid0] * (T(1) - midFrac) + cyL[iMid1] * midFrac;
    for (int i = 0; i < nSp; ++i) {
      cxL[i] -= xMid;
      cyL[i] -= yMid;
    }
  } else {
    for (int i = 0; i < nSp; ++i) {
      cxL[i] = scale * (s[i] - sCm);
      cyL[i] = scale * Aval[i] * std::sin(phase[i]);

      const T dxDs = scale;
      const T dyDs = scale * slope[i];
      const T tlen = std::sqrt(dxDs * dxDs + dyDs * dyDs);
      txV[i] = dxDs / tlen;
      tyV[i] = dyDs / tlen;
      nxN[i] = -tyV[i];
      nyN[i] =  txV[i];
    }
  }

  // Upper/lower surfaces in local frame
  std::vector<T> upperXL(nSp), upperYL(nSp), lowerXL(nSp), lowerYL(nSp);
  for (int i = 0; i < nSp; ++i) {
    upperXL[i] = cxL[i] + r * nxN[i];
    upperYL[i] = cyL[i] + r * nyN[i];
    lowerXL[i] = cxL[i] - r * nxN[i];
    lowerYL[i] = cyL[i] - r * nyN[i];
  }

  // Semicircular caps
  const int nHeadCap = std::max(static_cast<int>(std::ceil(M_PI * r / dsSpine)), 8);
  const int nTailCap = std::max(static_cast<int>(std::ceil(M_PI * r / dsSpine)), 8);
  const T dsHeadCap = M_PI * r / nHeadCap;
  const T dsTailCap = M_PI * r / nTailCap;

  // Head cap
  std::vector<T> headCapX(nHeadCap + 1), headCapY(nHeadCap + 1);
  for (int i = 0; i <= nHeadCap; ++i) {
    const T phi = M_PI * i / nHeadCap;
    headCapX[i] = cxL[0] + r * (-std::cos(phi) * nxN[0] - std::sin(phi) * txV[0]);
    headCapY[i] = cyL[0] + r * (-std::cos(phi) * nyN[0] - std::sin(phi) * tyV[0]);
  }
  // Tail cap
  std::vector<T> tailCapX(nTailCap + 1), tailCapY(nTailCap + 1);
  for (int i = 0; i <= nTailCap; ++i) {
    const T phi = M_PI * i / nTailCap;
    tailCapX[i] = cxL[nSp-1] + r * (std::cos(phi) * nxN[nSp-1] + std::sin(phi) * txV[nSp-1]);
    tailCapY[i] = cyL[nSp-1] + r * (std::cos(phi) * nyN[nSp-1] + std::sin(phi) * tyV[nSp-1]);
  }

  // Assemble closed contour (no duplicate junctions)
  // upper side + tail cap interior + lower side reversed + head cap interior
  const int nTotal = nSp + (nTailCap - 1) + nSp + (nHeadCap - 1);
  std::vector<T> locX(nTotal), locY(nTotal);
  allDs.resize(nTotal);

  int idx = 0;
  // Upper side
  for (int i = 0; i < nSp; ++i) {
    locX[idx] = upperXL[i]; locY[idx] = upperYL[i];
    allDs[idx] = dsSpine;
    idx++;
  }
  // Tail cap (interior points: 1..nTailCap-1)
  for (int i = 1; i < nTailCap; ++i) {
    locX[idx] = tailCapX[i]; locY[idx] = tailCapY[i];
    allDs[idx] = dsTailCap;
    idx++;
  }
  // Lower side (reversed)
  for (int i = nSp - 1; i >= 0; --i) {
    locX[idx] = lowerXL[i]; locY[idx] = lowerYL[i];
    allDs[idx] = dsSpine;
    idx++;
  }
  // Head cap (interior points: 1..nHeadCap-1)
  for (int i = 1; i < nHeadCap; ++i) {
    locX[idx] = headCapX[i]; locY[idx] = headCapY[i];
    allDs[idx] = dsHeadCap;
    idx++;
  }

  // Rotate to global frame
  const T cosT = std::cos(theta);
  const T sinT = std::sin(theta);

  globX.resize(nTotal);
  globY.resize(nTotal);

  for (int i = 0; i < nTotal; ++i) {
    globX[i] = cosT * locX[i] - sinT * locY[i] + xCm;
    globY[i] = sinT * locX[i] + cosT * locY[i] + yCm;
  }
}

}  // namespace

void buildCapsuleGeometry(
    const EelParams& p,
    T t, T xCm, T yCm, T theta,
    T dtLbm, T ampRamp,
    // outputs
    std::vector<T>& globX, std::vector<T>& globY,
    std::vector<T>& vxDefGlob, std::vector<T>& vyDefGlob,
    std::vector<T>& allDs)
{
  buildCapsulePositionFrame(p, t, xCm, yCm, theta, ampRamp,
                            globX, globY, allDs);

  const int nTotal = static_cast<int>(globX.size());
  vxDefGlob.assign(nTotal, T(0));
  vyDefGlob.assign(nTotal, T(0));
  if (!(dtLbm > T(0)) || nTotal == 0) {
    return;
  }

  // Desired IBM deformation velocity must be the time derivative of the
  // actual marker contour. A centered one-step finite difference keeps
  // surface and cap velocities consistent with tangent/normal motion.
  std::vector<T> xPlus, yPlus, dsPlus;
  std::vector<T> xMinus, yMinus, dsMinus;
  buildCapsulePositionFrame(p, t + dtLbm, xCm, yCm, theta, ampRamp,
                            xPlus, yPlus, dsPlus);
  buildCapsulePositionFrame(p, t - dtLbm, xCm, yCm, theta, ampRamp,
                            xMinus, yMinus, dsMinus);

  if (xPlus.size() != globX.size() || xMinus.size() != globX.size() ||
      yPlus.size() != globY.size() || yMinus.size() != globY.size()) {
    return;
  }

  for (int i = 0; i < nTotal; ++i) {
    vxDefGlob[i] = T(0.5) * (xPlus[i] - xMinus[i]);
    vyDefGlob[i] = T(0.5) * (yPlus[i] - yMinus[i]);
  }
}

std::array<T, 2> bodyForwardAxis(T theta)
{
  // theta rotates the head-to-tail spine axis; swimming forward points to the head.
  return {-std::cos(theta), -std::sin(theta)};
}

std::array<T, 2> bodyLateralAxis(T theta)
{
  const auto eForward = bodyForwardAxis(theta);
  return {-eForward[1], eForward[0]};
}
