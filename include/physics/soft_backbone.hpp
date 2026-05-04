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

SoftBackboneConfig makeSoftBackboneConfig(
    const EelParams& p,
    const PlanarRodSectionEstimate& rodSection,
    int nSegments);

SoftBackboneState makeStraightBackboneState(int nSegments, T theta0 = 0.0);

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
