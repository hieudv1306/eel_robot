#pragma once

#include "core/enums.hpp"
#include "core/params.hpp"
#include "core/types.hpp"

#include <string>

struct RunConfig {
  EelParams p;
  RunMode runMode = RunMode::Full;
  StudyMode studyMode = StudyMode::Standard;
  SimulationCase simCase = SimulationCase::SurgeOnly;
  WarmupMode warmupMode = WarmupMode::Rest;
  GaitNormalization gaitNormalization = GaitNormalization::Fixed;
  WallBoundary wallBoundary = WallBoundary::NoSlip;
  bool exportVelocity = true;
  bool exportVorticity = true;
  bool exportDiagnostics = true;
  bool exportBody = true;
  BodyKinematics bodyKinematics = BodyKinematics::PrescribedWave;
  bool softBackboneDynamics = false;
  SoftBackboneIntegrator softBackboneIntegrator = SoftBackboneIntegrator::Implicit;
  T softBackboneRelaxationTime = 0.05;
  T softBackboneFluidTorqueScale = 1.0;
  T softBackboneMaxAngleStep = 0.02;
  // Fraction of the slender-body theoretical added rotational inertia
  // lumped into segment inertia.  Defaults to 1.0 (full theoretical added
  // mass) which is required to keep partitioned-FSI dynamics stable when
  // the displaced fluid mass is comparable to the body mass.  Set to 0.0
  // to recover the structural-only inertia for diagnostic comparisons.
  T softBackboneAddedMassFrac = 1.0;
  bool softBackboneAbortOnInstability = true;
  T softBackboneAbortMeanSlip = 0.5;
  T softBackboneAbortMaxSlip = 50.0;
  int softBackboneAbortSaturatedFrames = 3;
  T alphaIBM = 1.0;
  int ibmIterations = 1;
  bool legacyKappaInputUsed = false;
  T legacyKappaInputValue = 0.0;
  T targetSt = -1.0;
  T referenceU = -1.0;
  T tailAmpRatioTarget = -1.0;
  T tCut = -1.0;
  std::string summaryCsv;
  std::string sensitivityCsv;
  std::string runTag;
  bool rawGeometryOverride = false;
  bool aspectGeometryOverride = false;
};

RunConfig parseCommandLine(int argc, char* argv[], const std::string& baseLogDir);
