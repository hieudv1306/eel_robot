#pragma once

#include "core/enums.hpp"
#include "core/params.hpp"
#include "core/types.hpp"
#include "physics/diagnostics.hpp"
#include "physics/markers.hpp"

#include <cstdint>
#include <string>
#include <vector>

struct CsvPathResolution {
  std::string path;
  bool headerMismatch = false;
  std::string fallbackPath;
};

struct CsvWriteResult {
  CsvPathResolution resolution;
  bool opened = false;
};

struct SummaryCsvInputs {
  EelParams p;
  SteadySummary steady;
  CycleConvergence cycleConv;
  SimulationCase simCase = SimulationCase::SurgeOnly;
  WarmupMode warmupMode = WarmupMode::Rest;
  GaitNormalization gaitNormalization = GaitNormalization::Fixed;
  WallBoundary wallBoundary = WallBoundary::NoSlip;
  bool initialPlacementClamped = false;
  std::uint64_t initialPlacementClampCount = 0;
  bool runtimeDomainClampHit = false;
  std::uint64_t runtimeDomainClampCount = 0;
  int nx = 0;
  int ny = 0;
  T tCut = 0;
  T bodyCenterlineLength = 0;
  T bodyLength = 0;
  T bodyWidth = 0;
  T mass = 0;
  T Ibody = 0;
  T effectiveEelFreq = 0;
  T effectiveEelA0 = 0;
  T alphaIBM = 0;
  bool legacyKappaInputUsed = false;
  int ibmIterations = 1;
};

struct VerificationCsvInputs {
  EelParams p;
  SteadySummary steady;
  CycleConvergence cycleConv;
  std::string runId;
  SimulationCase simCase = SimulationCase::SurgeOnly;
  StudyMode studyMode = StudyMode::Standard;
  WarmupMode warmupMode = WarmupMode::Rest;
  GaitNormalization gaitNormalization = GaitNormalization::Fixed;
  WallBoundary wallBoundary = WallBoundary::NoSlip;
  bool initialPlacementClamped = false;
  std::uint64_t initialPlacementClampCount = 0;
  bool speedClampHit = false;
  bool runtimeDomainClampHit = false;
  std::uint64_t speedClampCount = 0;
  std::uint64_t runtimeDomainClampCount = 0;
  std::uint64_t omegaClampCount = 0;
  int nx = 0;
  int ny = 0;
  T dtLbm = 0;
  T tCut = 0;
  T bodyCenterlineLength = 0;
  T bodyLength = 0;
  T bodyWidth = 0;
  T mass = 0;
  T Ibody = 0;
  T effectiveEelFreq = 0;
  T effectiveEelA0 = 0;
  T alphaIBM = 0;
  bool legacyKappaInputUsed = false;
  int ibmIterations = 1;
};

struct HistoryCsvData {
  const std::vector<T>* histT = nullptr;
  const std::vector<T>* histVx = nullptr;
  const std::vector<T>* histVy = nullptr;
  const std::vector<T>* histWz = nullptr;
  const std::vector<T>* histFx = nullptr;
  const std::vector<T>* histFy = nullptr;
  const std::vector<T>* histTz = nullptr;
  const std::vector<T>* histXcm = nullptr;
  const std::vector<T>* histYcm = nullptr;
  const std::vector<T>* histTheta = nullptr;
  const std::vector<T>* histPower = nullptr;
  const std::vector<T>* histEta = nullptr;
  const std::vector<T>* histUswim = nullptr;
  const std::vector<T>* histUlat = nullptr;
  const std::vector<T>* histForwardNetForce = nullptr;
  const std::vector<T>* histFlat = nullptr;
  const std::vector<T>* histRe = nullptr;
  const std::vector<T>* histSt = nullptr;
  const std::vector<T>* histMach = nullptr;
  const std::vector<T>* histUstar = nullptr;
  const std::vector<T>* histMeanSlip = nullptr;
  const std::vector<T>* histMaxSlip = nullptr;
  const std::vector<T>* histMeanMarkerForce = nullptr;
  const std::vector<T>* histMaxMarkerForce = nullptr;
  const std::vector<T>* histNormalizedSlip = nullptr;
  const std::vector<T>* histPowerRigid = nullptr;
  const std::vector<T>* histPowerDef = nullptr;
  const std::vector<T>* histMeanResidualSlip = nullptr;
  const std::vector<T>* histMaxResidualSlip = nullptr;
};

CsvPathResolution resolveCsvPath(const std::string& requestedPath,
                                 const std::string& expectedHeader,
                                 const std::string& schemaSuffix);
CsvWriteResult appendArSummaryCsv(const std::string& requestedPath,
                                   const SummaryCsvInputs& input);
CsvWriteResult appendVerificationSummaryCsv(const std::string& requestedPath,
                                             const VerificationCsvInputs& input);
std::string writeHistoryCsv(const std::string& runOutDir, const HistoryCsvData& history);
std::string writeFinalBodyCsv(const std::string& runOutDir, const LagrangianMarkers& markers);
