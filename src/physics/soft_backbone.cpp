#include "physics/soft_backbone.hpp"

#include "physics/gait.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

T wrapAngle(T angle)
{
  angle = std::fmod(angle + M_PI, T(2) * M_PI);
  if (angle < T(0)) angle += T(2) * M_PI;
  return angle - M_PI;
}

struct WaveSample {
  T yNorm = 0.0;
  T slope = 0.0;
};

WaveSample samplePreferredWave(const EelParams& p, T fraction, T t,
                               T externalAmpScale)
{
  const T L = std::max(p.eelLength, T(1e-12));
  const T s = std::clamp(fraction, T(0), T(1)) * L;
  const T lambda = std::max(p.eelLambda, T(1e-12));
  const auto act = actuationProfile(t, p.restTime, p.rampTime);
  const T ampScale = act.ampScale * externalAmpScale;
  const T waveSign = static_cast<T>(waveDirectionSign(p.waveDirection));
  const T phase =
    T(2) * M_PI * (waveSign * s / lambda - p.eelFreq * act.tActive);
  const T env = amplitudeEnvelope(s, L, p.eelA0) * ampScale;
  const T dEnvDs = p.eelA0 * T(6.6) * std::pow(s / L, T(2)) / L
                 * ampScale;

  WaveSample out;
  out.yNorm = env * std::sin(phase);
  out.slope = dEnvDs * std::sin(phase)
            + env * std::cos(phase) * (waveSign * T(2) * M_PI / lambda);
  return out;
}

std::vector<T> samplePreferredSegmentAngles(
    const EelParams& p,
    const SoftBackboneConfig& config,
    T t,
    T externalAmpScale)
{
  std::vector<T> theta(config.nSegments, T(0));
  if (!config.valid || config.nSegments <= 0) {
    return theta;
  }

  if (p.geometryKinematics == GeometryKinematics::InextensibleWave) {
    for (int i = 0; i < config.nSegments; ++i) {
      const T fraction = (T(i) + T(0.5)) / T(config.nSegments);
      theta[i] =
        std::atan(samplePreferredWave(p, fraction, t,
                                      externalAmpScale).slope);
    }
    return theta;
  }

  const int nNodes = config.nSegments + 1;
  std::vector<T> x(nNodes), y(nNodes);
  const T paramLength = std::max(p.eelLength, T(1e-12));
  const T physicalScale = config.centerlineLengthM / paramLength;
  for (int i = 0; i < nNodes; ++i) {
    const T fraction = T(i) / T(config.nSegments);
    const T s = fraction * paramLength;
    const WaveSample wave =
      samplePreferredWave(p, fraction, t, externalAmpScale);
    x[i] = physicalScale * (s - T(0.5) * paramLength);
    y[i] = physicalScale * wave.yNorm;
  }

  for (int i = 0; i < config.nSegments; ++i) {
    theta[i] = std::atan2(y[i + 1] - y[i], x[i + 1] - x[i]);
  }
  return theta;
}

SoftBackboneCurvatureProfile curvatureFromSegmentAngles(
    const SoftBackboneConfig& config,
    const std::vector<T>& theta)
{
  SoftBackboneCurvatureProfile out;
  if (!config.valid || config.nSegments < 2 ||
      static_cast<int>(theta.size()) != config.nSegments ||
      !(config.dsM > T(0))) {
    return out;
  }

  const int nJoints = config.nSegments - 1;
  out.sJointM.resize(nJoints);
  out.curvature.resize(nJoints);
  out.curvatureRate.assign(nJoints, T(0));
  for (int j = 0; j < nJoints; ++j) {
    out.sJointM[j] = (T(j) + T(1)) * config.dsM;
    out.curvature[j] = wrapAngle(theta[j + 1] - theta[j]) / config.dsM;
    out.maxAbsCurvature =
      std::max(out.maxAbsCurvature, std::abs(out.curvature[j]));
  }
  return out;
}

}  // namespace

SoftBackboneConfig makeSoftBackboneConfig(
    const EelParams& p,
    const PlanarRodSectionEstimate& rodSection,
    int nSegments)
{
  SoftBackboneConfig out;
  if (!rodSection.valid || nSegments < 2 ||
      !(rodSection.physicalLengthM > T(0)) ||
      !(p.totalGeometricLengthLU() > T(0)) ||
      !(p.centerlineLengthLU() > T(0)) ||
      !(rodSection.bendingStiffnessNm2 > T(0)) ||
      !(rodSection.massPerLengthKgM > T(0))) {
    return out;
  }

  out.valid = true;
  out.nSegments = nSegments;
  out.physicalLengthM = rodSection.physicalLengthM;
  out.centerlineLengthM =
    rodSection.physicalLengthM * p.centerlineLengthLU()
    / p.totalGeometricLengthLU();
  out.dsM = out.centerlineLengthM / T(nSegments);
  out.bendingStiffnessNm2 = rodSection.bendingStiffnessNm2;
  out.axialStiffnessN = rodSection.axialStiffnessN;
  out.massPerLengthKgM = rodSection.massPerLengthKgM;
  out.segmentMassKg = out.massPerLengthKgM * out.dsM;
  out.segmentRotationalInertiaKgM2 =
    out.segmentMassKg * out.dsM * out.dsM / T(12);
  out.jointAngleStiffnessNm = out.bendingStiffnessNm2 / out.dsM;
  out.dampingRatio = rodSection.dampingRatio;

  if (out.dampingRatio > T(0) &&
      out.jointAngleStiffnessNm > T(0) &&
      out.segmentRotationalInertiaKgM2 > T(0)) {
    const T effectiveRelativeInertia =
      T(0.5) * out.segmentRotationalInertiaKgM2;
    out.jointAngleDampingNms =
      T(2) * out.dampingRatio
      * std::sqrt(out.jointAngleStiffnessNm * effectiveRelativeInertia);
    out.curvatureDampingNm2s = out.jointAngleDampingNms * out.dsM;
  }
  return out;
}

SoftBackboneState makeStraightBackboneState(int nSegments, T theta0)
{
  SoftBackboneState state;
  if (nSegments <= 0) {
    return state;
  }
  state.theta.assign(nSegments, theta0);
  state.omega.assign(nSegments, T(0));
  return state;
}

SoftBackboneState preferredBackboneStateWave(
    const EelParams& p,
    const SoftBackboneConfig& config,
    T t,
    T rateDt,
    T externalAmpScale)
{
  SoftBackboneState state;
  if (!config.valid || config.nSegments <= 0) {
    return state;
  }

  state.theta =
    samplePreferredSegmentAngles(p, config, t, externalAmpScale);
  state.omega.assign(state.theta.size(), T(0));
  if (!(rateDt > T(0)) || state.theta.empty()) {
    return state;
  }

  const std::vector<T> plus =
    samplePreferredSegmentAngles(p, config, t + rateDt, externalAmpScale);
  const std::vector<T> minus =
    samplePreferredSegmentAngles(p, config, t - rateDt, externalAmpScale);
  if (plus.size() != state.theta.size() ||
      minus.size() != state.theta.size()) {
    return state;
  }

  for (size_t i = 0; i < state.theta.size(); ++i) {
    state.omega[i] = wrapAngle(plus[i] - minus[i]) / (T(2) * rateDt);
  }
  return state;
}

SoftBackboneCurvatureProfile preferredBackboneCurvatureWave(
    const EelParams& p,
    const SoftBackboneConfig& config,
    T t,
    T rateDt,
    T externalAmpScale)
{
  SoftBackboneCurvatureProfile out =
    curvatureFromSegmentAngles(
      config, samplePreferredSegmentAngles(p, config, t, externalAmpScale));
  if (!config.valid || !(rateDt > T(0)) || out.curvature.empty()) {
    return out;
  }

  SoftBackboneCurvatureProfile plus =
    curvatureFromSegmentAngles(config,
                               samplePreferredSegmentAngles(
                                 p, config, t + rateDt, externalAmpScale));
  SoftBackboneCurvatureProfile minus =
    curvatureFromSegmentAngles(config,
                               samplePreferredSegmentAngles(
                                 p, config, t - rateDt, externalAmpScale));
  if (plus.curvature.size() != out.curvature.size() ||
      minus.curvature.size() != out.curvature.size()) {
    return out;
  }

  out.curvatureRate.assign(out.curvature.size(), T(0));
  for (size_t i = 0; i < out.curvature.size(); ++i) {
    out.curvatureRate[i] =
      (plus.curvature[i] - minus.curvature[i]) / (T(2) * rateDt);
    out.maxAbsCurvatureRate =
      std::max(out.maxAbsCurvatureRate, std::abs(out.curvatureRate[i]));
  }
  return out;
}

SoftBackboneCurvatureProfile curvatureFromBackboneState(
    const SoftBackboneConfig& config,
    const SoftBackboneState& state)
{
  SoftBackboneCurvatureProfile out =
    curvatureFromSegmentAngles(config, state.theta);
  if (!config.valid ||
      state.omega.size() != state.theta.size() ||
      out.curvature.empty() ||
      !(config.dsM > T(0))) {
    return out;
  }

  out.curvatureRate.assign(out.curvature.size(), T(0));
  for (size_t j = 0; j < out.curvature.size(); ++j) {
    out.curvatureRate[j] =
      (state.omega[j + 1] - state.omega[j]) / config.dsM;
    out.maxAbsCurvatureRate =
      std::max(out.maxAbsCurvatureRate, std::abs(out.curvatureRate[j]));
  }
  return out;
}

SoftBackboneTorqueResult computeSoftBackboneTorques(
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    const SoftBackboneCurvatureProfile& preferred)
{
  SoftBackboneTorqueResult out;
  const int nSegments = config.nSegments;
  const int nJoints = nSegments - 1;
  if (!config.valid || nSegments < 2 ||
      static_cast<int>(state.theta.size()) != nSegments ||
      static_cast<int>(state.omega.size()) != nSegments ||
      static_cast<int>(preferred.curvature.size()) != nJoints ||
      !(config.dsM > T(0))) {
    return out;
  }

  out.jointMomentNm.assign(nJoints, T(0));
  out.segmentTorqueNm.assign(nSegments, T(0));
  for (int j = 0; j < nJoints; ++j) {
    const T actualCurvature =
      wrapAngle(state.theta[j + 1] - state.theta[j]) / config.dsM;
    const T actualCurvatureRate =
      (state.omega[j + 1] - state.omega[j]) / config.dsM;
    const T preferredCurvature = preferred.curvature[j];
    const T preferredRate =
      (j < static_cast<int>(preferred.curvatureRate.size()))
        ? preferred.curvatureRate[j] : T(0);

    const T curvatureError = actualCurvature - preferredCurvature;
    const T rateError = actualCurvatureRate - preferredRate;
    const T moment =
      config.bendingStiffnessNm2 * curvatureError
      + config.curvatureDampingNm2s * rateError;

    out.jointMomentNm[j] = moment;
    out.segmentTorqueNm[j] += moment;
    out.segmentTorqueNm[j + 1] -= moment;
    out.elasticEnergyJ +=
      T(0.5) * config.bendingStiffnessNm2 * config.dsM
      * curvatureError * curvatureError;
    out.maxAbsJointMomentNm =
      std::max(out.maxAbsJointMomentNm, std::abs(moment));
  }

  for (T torque : out.segmentTorqueNm) {
    out.maxAbsSegmentTorqueNm =
      std::max(out.maxAbsSegmentTorqueNm, std::abs(torque));
  }
  return out;
}

SoftBackboneCenterline buildSoftBackboneCenterlineLU(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    T xCm,
    T yCm,
    T theta)
{
  SoftBackboneCenterline out;
  if (!config.valid || config.nSegments <= 0 ||
      static_cast<int>(state.theta.size()) != config.nSegments ||
      !(config.centerlineLengthM > T(0))) {
    return out;
  }

  const T scaleToLU = p.centerlineLengthLU() / config.centerlineLengthM;
  const T dsLu = config.dsM * scaleToLU;
  if (!(dsLu > T(0))) {
    return out;
  }

  const int nSegments = config.nSegments;
  const int nNodes = nSegments + 1;
  std::vector<T> localX(nNodes, T(0));
  std::vector<T> localY(nNodes, T(0));
  for (int i = 0; i < nSegments; ++i) {
    localX[i + 1] = localX[i] + dsLu * std::cos(state.theta[i]);
    localY[i + 1] = localY[i] + dsLu * std::sin(state.theta[i]);
  }

  const T midSeg = T(0.5) * T(nSegments);
  const int iMid0 = static_cast<int>(std::floor(midSeg));
  const int iMid1 = std::min(iMid0 + 1, nNodes - 1);
  const T midFrac = midSeg - T(iMid0);
  const T xMid = localX[iMid0] * (T(1) - midFrac) + localX[iMid1] * midFrac;
  const T yMid = localY[iMid0] * (T(1) - midFrac) + localY[iMid1] * midFrac;

  const T cosT = std::cos(theta);
  const T sinT = std::sin(theta);
  out.nodeX.resize(nNodes);
  out.nodeY.resize(nNodes);
  for (int i = 0; i < nNodes; ++i) {
    const T xLocal = localX[i] - xMid;
    const T yLocal = localY[i] - yMid;
    out.nodeX[i] = cosT * xLocal - sinT * yLocal + xCm;
    out.nodeY[i] = sinT * xLocal + cosT * yLocal + yCm;
  }

  out.segmentX.resize(nSegments);
  out.segmentY.resize(nSegments);
  out.segmentTheta.resize(nSegments);
  for (int i = 0; i < nSegments; ++i) {
    out.segmentX[i] = T(0.5) * (out.nodeX[i] + out.nodeX[i + 1]);
    out.segmentY[i] = T(0.5) * (out.nodeY[i] + out.nodeY[i + 1]);
    out.segmentTheta[i] = state.theta[i] + theta;
  }
  return out;
}

SoftBackboneState extrapolateBackboneState(
    const SoftBackboneState& state,
    T dt)
{
  SoftBackboneState out = state;
  if (out.theta.size() != out.omega.size()) {
    return out;
  }
  for (size_t i = 0; i < out.theta.size(); ++i) {
    out.theta[i] = wrapAngle(out.theta[i] + out.omega[i] * dt);
  }
  return out;
}

SoftBackboneDynamicsDiagnostics advanceSoftBackboneOverdamped(
    const SoftBackboneConfig& config,
    const SoftBackboneState& preferred,
    const std::vector<T>& fluidSegmentTorqueNm,
    T dt,
    const SoftBackboneDynamicsParams& params,
    SoftBackboneState& state)
{
  SoftBackboneDynamicsDiagnostics diag;
  const int nSegments = config.nSegments;
  if (!config.valid || nSegments < 2 || !(dt > T(0)) ||
      static_cast<int>(state.theta.size()) != nSegments ||
      static_cast<int>(state.omega.size()) != nSegments ||
      static_cast<int>(preferred.theta.size()) != nSegments ||
      static_cast<int>(preferred.omega.size()) != nSegments ||
      static_cast<int>(fluidSegmentTorqueNm.size()) != nSegments ||
      !(config.jointAngleStiffnessNm > T(0))) {
    return diag;
  }

  const T relax = std::max(params.relaxationTime, dt);
  const T blend = std::clamp(dt / relax, T(0), T(1));
  const T maxStep = (params.maxAngleStep > T(0))
                  ? params.maxAngleStep : T(1e30);

  std::vector<T> targetJointAngle(nSegments - 1, T(0));
  for (int j = 0; j < nSegments - 1; ++j) {
    const T qPreferred = wrapAngle(preferred.theta[j + 1] -
                                   preferred.theta[j]);
    const T fluidGeneralizedTorque =
      params.fluidTorqueScale *
      (fluidSegmentTorqueNm[j + 1] - fluidSegmentTorqueNm[j]);
    diag.maxAbsFluidSegmentTorqueNm =
      std::max(diag.maxAbsFluidSegmentTorqueNm,
               std::max(std::abs(fluidSegmentTorqueNm[j]),
                        std::abs(fluidSegmentTorqueNm[j + 1])));
    const T curvatureOffset =
      fluidGeneralizedTorque / config.bendingStiffnessNm2;
    diag.maxAbsTargetCurvatureOffset =
      std::max(diag.maxAbsTargetCurvatureOffset,
               std::abs(curvatureOffset));
    targetJointAngle[j] =
      qPreferred + curvatureOffset * config.dsM;
  }

  // Pin segment 0 to the preferred state.  This removes the rigid-rotation
  // null space (uniform offset in theta is a kinematic invariant of the
  // joint-angle stiffness operator) without the meanError drift correction
  // that the previous version applied to the whole segment array.  The
  // rigid-body theta absorbs net rotation separately.
  std::vector<T> targetTheta(nSegments, T(0));
  targetTheta[0] = preferred.theta[0];
  for (int i = 1; i < nSegments; ++i) {
    targetTheta[i] = targetTheta[i - 1] + targetJointAngle[i - 1];
  }

  state.theta[0] = preferred.theta[0];
  state.omega[0] = preferred.omega[0];
  diag.maxAbsAngleStep = std::max(diag.maxAbsAngleStep,
                                   std::abs(preferred.omega[0] * dt));

  for (int i = 1; i < nSegments; ++i) {
    T step = blend * wrapAngle(targetTheta[i] - state.theta[i]);
    step = std::clamp(step, -maxStep, maxStep);
    state.omega[i] = step / dt;
    state.theta[i] = wrapAngle(state.theta[i] + step);
    diag.maxAbsAngleStep = std::max(diag.maxAbsAngleStep, std::abs(step));
    diag.maxAbsAngleErrorRad =
      std::max(diag.maxAbsAngleErrorRad,
               std::abs(wrapAngle(state.theta[i] - preferred.theta[i])));
  }
  return diag;
}

SoftBackboneDynamicsDiagnostics advanceSoftBackboneImplicit(
    const SoftBackboneConfig& config,
    const SoftBackboneState& preferred,
    const std::vector<T>& fluidSegmentTorqueNm,
    T dt,
    const SoftBackboneDynamicsParams& params,
    SoftBackboneState& state)
{
  SoftBackboneDynamicsDiagnostics diag;
  const int nSegments = config.nSegments;
  if (!config.valid || nSegments < 2 || !(dt > T(0)) ||
      static_cast<int>(state.theta.size()) != nSegments ||
      static_cast<int>(state.omega.size()) != nSegments ||
      static_cast<int>(preferred.theta.size()) != nSegments ||
      static_cast<int>(preferred.omega.size()) != nSegments ||
      static_cast<int>(fluidSegmentTorqueNm.size()) != nSegments ||
      !(config.jointAngleStiffnessNm > T(0)) ||
      !(config.segmentRotationalInertiaKgM2 > T(0))) {
    return diag;
  }

  const int nInterior = nSegments - 1;  // segments 1..N-1
  const int nJoints = nSegments - 1;    // joints 0..N-2

  const T I = config.segmentRotationalInertiaKgM2;
  const T K = config.jointAngleStiffnessNm;
  const T C = config.jointAngleDampingNms;
  const T A = K * dt + C;
  const T diagBase = I / dt;
  const T fluidScale = params.fluidTorqueScale;
  const T maxStep = (params.maxAngleStep > T(0))
                  ? params.maxAngleStep : T(1e30);

  // Track raw fluid torque magnitudes for diagnostics (independent of scale).
  for (int i = 0; i < nSegments; ++i) {
    diag.maxAbsFluidSegmentTorqueNm =
      std::max(diag.maxAbsFluidSegmentTorqueNm,
               std::abs(fluidSegmentTorqueNm[i]));
  }

  // RHS_j = K (theta_{j+1}^n - theta_j^n - q_pref_j) - C * qdot_pref_j,
  // i.e. the part of joint moment M_j^{n+1} that is independent of omega^{n+1}
  // after substituting theta^{n+1} = theta^n + dt * omega^{n+1}.
  std::vector<T> rhsJoint(nJoints, T(0));
  for (int j = 0; j < nJoints; ++j) {
    const T qPref = wrapAngle(preferred.theta[j + 1] - preferred.theta[j]);
    const T qDotPref = preferred.omega[j + 1] - preferred.omega[j];
    const T thetaDiff = wrapAngle(state.theta[j + 1] - state.theta[j]);
    rhsJoint[j] = K * (thetaDiff - qPref) - C * qDotPref;
    diag.maxAbsTargetCurvatureOffset =
      std::max(diag.maxAbsTargetCurvatureOffset,
               std::abs(qPref) / std::max(config.dsM, T(1e-30)));
  }

  // Pin segment 0 kinematically.
  const T omega0Pinned = preferred.omega[0];
  const T theta0Pinned = preferred.theta[0];

  // Solve tridiagonal system for omega_i^{n+1}, i = 1..N-1.
  // Variables indexed locally as k = i - 1, k = 0..nInterior-1.
  //   sub_k * omega_{k-1} + diag_k * omega_k + sup_k * omega_{k+1} = rhs_k
  std::vector<T> sub(nInterior, T(0));
  std::vector<T> dia(nInterior, T(0));
  std::vector<T> sup(nInterior, T(0));
  std::vector<T> rhs(nInterior, T(0));

  for (int k = 0; k < nInterior; ++k) {
    const int i = k + 1;
    const bool isLast = (i == nSegments - 1);
    const T fluid = fluidScale * fluidSegmentTorqueNm[i];

    if (isLast) {
      // Segment N-1: tau_int = -M_{N-2}, only.
      //   I/dt omega_{N-1} + A (omega_{N-1} - omega_{N-2}) = I/dt omega_{N-1}^n
      //                                                    + tau_ext - rhsJoint_{N-2}
      sub[k] = -A;
      dia[k] = diagBase + A;
      sup[k] = T(0);
      rhs[k] = diagBase * state.omega[i] + fluid - rhsJoint[i - 1];
    } else {
      // Interior segment: tau_int = M_i - M_{i-1}.
      //   -A omega_{i-1} + (I/dt + 2A) omega_i - A omega_{i+1}
      //     = I/dt omega_i^n + tau_ext + rhsJoint_i - rhsJoint_{i-1}
      sub[k] = -A;
      dia[k] = diagBase + T(2) * A;
      sup[k] = -A;
      rhs[k] = diagBase * state.omega[i] + fluid
             + rhsJoint[i] - rhsJoint[i - 1];
    }
  }

  // Apply the pinned-omega boundary at k=0 (i=1): omega_0 is known, so the
  // -A * omega_0 contribution moves to the RHS.
  if (nInterior > 0) {
    rhs[0] += A * omega0Pinned;
    sub[0] = T(0);
  }

  // Thomas elimination (in-place on dia/sup/rhs).
  for (int k = 1; k < nInterior; ++k) {
    const T denom = dia[k - 1];
    if (!(std::abs(denom) > T(1e-30))) {
      return diag;  // singular — leave state untouched
    }
    const T m = sub[k] / denom;
    dia[k] -= m * sup[k - 1];
    rhs[k] -= m * rhs[k - 1];
  }
  std::vector<T> omegaNew(nInterior, T(0));
  if (nInterior > 0) {
    if (!(std::abs(dia[nInterior - 1]) > T(1e-30))) {
      return diag;
    }
    omegaNew[nInterior - 1] = rhs[nInterior - 1] / dia[nInterior - 1];
    for (int k = nInterior - 2; k >= 0; --k) {
      if (!(std::abs(dia[k]) > T(1e-30))) {
        return diag;
      }
      omegaNew[k] = (rhs[k] - sup[k] * omegaNew[k + 1]) / dia[k];
    }
  }

  // Apply the solution.  Clamp |dt * omega| to maxStep as a numerical
  // safety; this should normally be inactive once the coupling is healthy.
  state.theta[0] = theta0Pinned;
  state.omega[0] = omega0Pinned;
  diag.maxAbsAngleStep = std::max(diag.maxAbsAngleStep,
                                   std::abs(omega0Pinned * dt));

  for (int k = 0; k < nInterior; ++k) {
    const int i = k + 1;
    T step = omegaNew[k] * dt;
    if (std::abs(step) > maxStep) {
      step = (step > T(0)) ? maxStep : -maxStep;
      omegaNew[k] = step / dt;
    }
    state.omega[i] = omegaNew[k];
    state.theta[i] = wrapAngle(state.theta[i] + step);
    diag.maxAbsAngleStep = std::max(diag.maxAbsAngleStep, std::abs(step));
    diag.maxAbsAngleErrorRad =
      std::max(diag.maxAbsAngleErrorRad,
               std::abs(wrapAngle(state.theta[i] - preferred.theta[i])));
  }
  return diag;
}
