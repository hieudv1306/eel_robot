#pragma once

#include "core/enums.hpp"
#include "core/params.hpp"
#include "core/types.hpp"

#include <string>
#include <vector>

struct RunConfig {
  EelParams p;
  RunMode runMode = RunMode::Full;
  StudyMode studyMode = StudyMode::Standard;
  SimulationCase simCase = SimulationCase::SurgeOnly;
  WarmupMode warmupMode = WarmupMode::Rest;
  WallBoundary wallBoundary = WallBoundary::NoSlip;
  bool exportVelocity = false;
  bool exportVorticity = false;
  bool exportDiagnostics = false;
  bool exportBody = false;
  BodyKinematics bodyKinematics = BodyKinematics::SoftBackbone;
  bool softBackboneDynamics = false;
  T softBackboneRelaxationTime = 0.05;
  T softBackboneFluidTorqueScale = 1.0;
  T softBackboneFluidTorqueFilterTime = 0.0;
  T softBackboneMaxAngleStep = 0.02;
  int softBackboneCouplingIterations = 1;
  T softBackboneCouplingRelaxation = 1.0;
  T softBackboneCouplingTolerance = 1e-4;
  SoftBackboneLoadProjection softBackboneLoadProjection =
    SoftBackboneLoadProjection::SegmentCentroid;
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
  int legacyNIbmIters = 1;
  T tCut = -1.0;
  std::string summaryCsv;
  std::string sensitivityCsv;
  std::string runTag;
  bool rawGeometryOverride = false;
  bool aspectGeometryOverride = false;
};

struct CliDiagnostics {
  std::vector<std::string> errors;
  std::vector<std::string> warnings;

  bool ok() const { return errors.empty(); }
};

struct CliParseResult {
  RunConfig config;
  CliDiagnostics diagnostics;
  bool helpRequested = false;
};

CliParseResult parseCommandLineDetailed(int argc, char* argv[],
                                        const std::string& baseLogDir);
RunConfig parseCommandLine(int argc, char* argv[], const std::string& baseLogDir);
CliDiagnostics validateRunConfig(const RunConfig& config);
std::string commandLineUsage(const char* executableName);
