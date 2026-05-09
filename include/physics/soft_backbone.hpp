#pragma once

#include "core/params.hpp"
#include "core/types.hpp"
#include "physics/soft_rod.hpp"

#include <vector>

struct SoftBackboneConfig {
  bool valid = false;
  int nSegments = 0;
  T physicalLengthM = 0.0;
  T centerlineLengthM = 0.0;
  T dsM = 0.0;
  T bendingStiffnessNm2 = 0.0;
  T axialStiffnessN = 0.0;
  T massPerLengthKgM = 0.0;
  T segmentMassKg = 0.0;
  T segmentRotationalInertiaKgM2 = 0.0;
  T jointAngleStiffnessNm = 0.0;
  T jointAngleDampingNms = 0.0;
  T curvatureDampingNm2s = 0.0;
  T dampingRatio = 0.0;
};

struct SoftBackboneState {
  std::vector<T> theta;
  std::vector<T> omega;
};

struct SoftBackboneCurvatureProfile {
  std::vector<T> sJointM;
  std::vector<T> curvature;
  std::vector<T> curvatureRate;
  T maxAbsCurvature = 0.0;
  T maxAbsCurvatureRate = 0.0;
};

struct SoftBackboneTorqueResult {
  std::vector<T> jointMomentNm;
  std::vector<T> segmentTorqueNm;
  T elasticEnergyJ = 0.0;
  T maxAbsJointMomentNm = 0.0;
  T maxAbsSegmentTorqueNm = 0.0;
};

struct SoftBackboneCenterline {
  std::vector<T> nodeX;
  std::vector<T> nodeY;
  std::vector<T> segmentX;
  std::vector<T> segmentY;
  std::vector<T> segmentTheta;
};

struct SoftBackboneDynamicsParams {
  T relaxationTime = 0.05;
  T fluidTorqueScale = 1.0;
  T maxAngleStep = 0.02;
};

struct SoftBackboneDynamicsDiagnostics {
  T maxAbsFluidSegmentTorqueNm = 0.0;
  T maxAbsTargetCurvatureOffset = 0.0;
  T maxAbsAngleStep = 0.0;
  T maxAbsAngleErrorRad = 0.0;
};

SoftBackboneConfig makeSoftBackboneConfig(
    const EelParams& p,
    const PlanarRodSectionEstimate& rodSection,
    int nSegments);

SoftBackboneState makeStraightBackboneState(int nSegments, T theta0 = 0.0);

// Derive nNodes = nSegments+1 node tangent angles from nSegments segment
// angles using centered averaging at interior nodes and linear extrapolation
// at the head/tail.  This makes the segment-angle representation produce a
// node-tangent sampling that is compatible with the legacy node-based
// inextensible-wave centerline (which samples slope analytically at every
// node and walks the centerline with mid-of-node-tangent / trapezoidal rule).
// Required for the soft-backbone marker geometry to match the legacy
// prescribed-wave geometry to O(ds^2) at the body endpoints, where the slope
// changes most rapidly.
std::vector<T> nodeTangentsFromSegmentAngles(
    const std::vector<T>& segmentAngles);

SoftBackboneState preferredBackboneStateWave(
    const EelParams& p,
    const SoftBackboneConfig& config,
    T t,
    T rateDt = 0.0,
    T externalAmpScale = 1.0);

SoftBackboneCurvatureProfile preferredBackboneCurvatureWave(
    const EelParams& p,
    const SoftBackboneConfig& config,
    T t,
    T rateDt = 0.0,
    T externalAmpScale = 1.0);

SoftBackboneCurvatureProfile curvatureFromBackboneState(
    const SoftBackboneConfig& config,
    const SoftBackboneState& state);

SoftBackboneTorqueResult computeSoftBackboneTorques(
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    const SoftBackboneCurvatureProfile& preferred);

SoftBackboneCenterline buildSoftBackboneCenterlineLU(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    T xCm,
    T yCm,
    T theta);

SoftBackboneState extrapolateBackboneState(
    const SoftBackboneState& state,
    T dt);

SoftBackboneDynamicsDiagnostics advanceSoftBackboneOverdamped(
    const SoftBackboneConfig& config,
    const SoftBackboneState& preferred,
    const std::vector<T>& fluidSegmentTorqueNm,
    T dt,
    const SoftBackboneDynamicsParams& params,
    SoftBackboneState& state);

// Implicit Euler advance of segment angles with proper inertia, joint
// stiffness K_theta = EI/ds, and joint damping C_theta from the material
// damping ratio.  Pins segment 0 to the preferred state to remove the
// rigid-rotation null space (the rigid-body theta absorbs net rotation).
// The system is tridiagonal in (N-1) interior segment angular velocities
// and is solved with Thomas elimination, so cost is O(N) per step and the
// scheme is unconditionally stable for stiff K_theta.
SoftBackboneDynamicsDiagnostics advanceSoftBackboneImplicit(
    const SoftBackboneConfig& config,
    const SoftBackboneState& preferred,
    const std::vector<T>& fluidSegmentTorqueNm,
    T dt,
    const SoftBackboneDynamicsParams& params,
    SoftBackboneState& state);
