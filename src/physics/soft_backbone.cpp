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
  const T phase = T(2) * M_PI * (s / lambda - p.eelFreq * act.tActive);
  const T env = amplitudeEnvelope(s, L, p.eelA0) * ampScale;
  const T dEnvDs = p.eelA0 * T(6.6) * std::pow(s / L, T(2)) / L
                 * ampScale;

  WaveSample out;
  out.yNorm = env * std::sin(phase);
  out.slope = dEnvDs * std::sin(phase)
            + env * std::cos(phase) * (T(2) * M_PI / lambda);
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
