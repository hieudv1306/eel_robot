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
  // Slender-body added rotational inertia per segment, scaled by
  // softBackboneAddedMassFrac.  Lumped into the structural inertia in
  // the implicit integrator to suppress the partitioned-FSI added-mass
  // instability (Causin et al. 2005) that otherwise appears when the
  // displaced fluid mass is comparable to the body mass.
  T addedSegmentRotationalInertiaKgM2 = 0.0;
  T effectiveSegmentRotationalInertiaKgM2 = 0.0;
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
  T fluidTorqueScale = 1.0;
  T maxAngleStep = 0.02;
};

struct SoftBackboneDynamicsDiagnostics {
  T maxAbsFluidSegmentTorqueNm = 0.0;
  T maxAbsTargetCurvatureOffset = 0.0;
  T maxAbsAngleStep = 0.0;
  T maxAbsAngleErrorRad = 0.0;
  T elasticEnergyJ = 0.0;
  T dampingPowerW = 0.0;
  // Signed proxy for work needed to move the preferred/rest joint angles
  // against the current spring-damper load.  Positive values mean the
  // prescribed soft backbone is doing work on the elastic model under this
  // sign convention; use absActuatorPowerProxyW for ranking.
  T actuatorPowerProxyW = 0.0;
  T absActuatorPowerProxyW = 0.0;
  // Power from the scaled fluid torque that is actually applied to the
  // backbone dynamics.  The raw fluid torque magnitude is still reported by
  // maxAbsFluidSegmentTorqueNm.
  T appliedFluidPowerW = 0.0;
  T absAppliedFluidPowerW = 0.0;
  T maxAbsJointMomentNm = 0.0;
};

SoftBackboneConfig makeSoftBackboneConfig(
    const EelParams& p,
    const PlanarRodSectionEstimate& rodSection,
    int nSegments,
    T addedMassFrac = T(1));

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

// Implicit Euler advance of segment angles with proper inertia, joint
// stiffness K_theta = EI/ds, and joint damping C_theta from the material
// damping ratio.  The rod is solved with free-end internal torques and the
// rigid-rotation mode is removed by matching the mean angle/rate to the
// preferred state (the rigid-body theta absorbs net rotation).  The system is
// tridiagonal in N segment angular velocities and is solved with Thomas
// elimination, so cost is O(N) per step and the scheme is unconditionally
// stable for stiff K_theta.
SoftBackboneDynamicsDiagnostics advanceSoftBackboneImplicit(
    const SoftBackboneConfig& config,
    const SoftBackboneState& preferred,
    const std::vector<T>& fluidSegmentTorqueNm,
    T dt,
    const SoftBackboneDynamicsParams& params,
    SoftBackboneState& state);
