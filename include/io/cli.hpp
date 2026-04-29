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
