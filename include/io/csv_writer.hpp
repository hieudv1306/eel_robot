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
  WallBoundary wallBoundary = WallBoundary::NoSlip;
  BodyKinematics bodyKinematics = BodyKinematics::SoftBackbone;
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
  T alphaIBM = 0;
  int ibmIterations = 1;
  bool softBackboneDynamics = false;
  T softBackboneRelaxationTime = 0;
  T softBackboneFluidTorqueScale = 0;
  T softBackboneFluidTorqueFilterTime = 0;
  T softBackboneAddedMassFrac = 0;
  T softBackboneMaxAngleStep = 0;
  int softBackboneCouplingIterations = 1;
  T softBackboneCouplingRelaxation = 1;
  T softBackboneCouplingTolerance = 0;
  SoftBackboneLoadProjection softBackboneLoadProjection =
    SoftBackboneLoadProjection::SegmentCentroid;
  T meanSoftFluidTorqueNm = 0;
  T maxSoftFluidTorqueNm = 0;
  T meanSoftAngleStep = 0;
  T maxSoftAngleStep = 0;
  T meanSoftCouplingResidual = 0;
  T maxSoftCouplingResidual = 0;
  T meanSoftCouplingItersUsed = 0;
  T maxSoftCouplingItersUsed = 0;
  T meanSoftElasticEnergyJ = 0;
  T maxSoftElasticEnergyJ = 0;
  T meanSoftDampingPowerW = 0;
  T meanSoftActuatorPowerProxyW = 0;
  T meanAbsSoftActuatorPowerProxyW = 0;
  T maxAbsSoftActuatorPowerProxyW = 0;
  T meanSoftAppliedFluidPowerW = 0;
  T meanAbsSoftAppliedFluidPowerW = 0;
  T maxAbsSoftAppliedFluidPowerW = 0;
  T softBodyMassKg = 0;
  T meanUPhysicalMps = 0;
  T cotSoftActuatorProxySI = 0;
};

struct VerificationCsvInputs {
  EelParams p;
  SteadySummary steady;
  CycleConvergence cycleConv;
  std::string runId;
  SimulationCase simCase = SimulationCase::SurgeOnly;
  StudyMode studyMode = StudyMode::Standard;
  WarmupMode warmupMode = WarmupMode::Rest;
  WallBoundary wallBoundary = WallBoundary::NoSlip;
  BodyKinematics bodyKinematics = BodyKinematics::SoftBackbone;
  bool initialPlacementClamped = false;
  std::uint64_t initialPlacementClampCount = 0;
  bool runtimeDomainClampHit = false;
  std::uint64_t runtimeDomainClampCount = 0;
  int nx = 0;
  int ny = 0;
  T dtLbm = 0;
  T tCut = 0;
  T bodyCenterlineLength = 0;
  T bodyLength = 0;
  T bodyWidth = 0;
  T mass = 0;
  T Ibody = 0;
  T alphaIBM = 0;
  int ibmIterations = 1;
  bool softBackboneDynamics = false;
  T softBackboneRelaxationTime = 0;
  T softBackboneFluidTorqueScale = 0;
  T softBackboneFluidTorqueFilterTime = 0;
  T softBackboneAddedMassFrac = 0;
  T softBackboneMaxAngleStep = 0;
  int softBackboneCouplingIterations = 1;
  T softBackboneCouplingRelaxation = 1;
  T softBackboneCouplingTolerance = 0;
  SoftBackboneLoadProjection softBackboneLoadProjection =
    SoftBackboneLoadProjection::SegmentCentroid;
  T meanSoftFluidTorqueNm = 0;
  T maxSoftFluidTorqueNm = 0;
  T meanSoftAngleStep = 0;
  T maxSoftAngleStep = 0;
  T meanSoftCouplingResidual = 0;
  T maxSoftCouplingResidual = 0;
  T meanSoftCouplingItersUsed = 0;
  T maxSoftCouplingItersUsed = 0;
  T meanSoftElasticEnergyJ = 0;
  T maxSoftElasticEnergyJ = 0;
  T meanSoftDampingPowerW = 0;
  T meanSoftActuatorPowerProxyW = 0;
  T meanAbsSoftActuatorPowerProxyW = 0;
  T maxAbsSoftActuatorPowerProxyW = 0;
  T meanSoftAppliedFluidPowerW = 0;
  T meanAbsSoftAppliedFluidPowerW = 0;
  T maxAbsSoftAppliedFluidPowerW = 0;
  T softBodyMassKg = 0;
  T meanUPhysicalMps = 0;
  T cotSoftActuatorProxySI = 0;
};

struct HistoryFrameSample {
  T t = 0;
  T Vx = 0;
  T Vy = 0;
  T omegaZ = 0;
  T Fx = 0;
  T Fy = 0;
  T Tz = 0;
  T xCm = 0;
  T yCm = 0;
  T theta = 0;
  T power = 0;
  T eta = 0;
  T Uswim = 0;
  T Ulateral = 0;
  T forwardNetForce = 0;
  T Flat = 0;
  T Re = 0;
  T St = 0;
  T Mach = 0;
  T Ustar = 0;
  T meanSlip = 0;
  T maxSlip = 0;
  T meanMarkerForce = 0;
  T maxMarkerForce = 0;
  T normalizedSlip = 0;
  T powerRigid = 0;
  T powerDef = 0;
  T meanResidualSlip = 0;
  T maxResidualSlip = 0;
  T meanSoftFluidTorqueNm = 0;
  T maxSoftFluidTorqueNm = 0;
  T meanSoftAngleStep = 0;
  T maxSoftAngleStep = 0;
  T softCouplingResidual = 0;
  T maxSoftCouplingResidual = 0;
  T softCouplingItersUsed = 0;
  T maxSoftCouplingItersUsed = 0;
  T softElasticEnergyJ = 0;
  T maxSoftElasticEnergyJ = 0;
  T softDampingPowerW = 0;
  T softActuatorPowerProxyW = 0;
  T absSoftActuatorPowerProxyW = 0;
  T maxAbsSoftActuatorPowerProxyW = 0;
  T softAppliedFluidPowerW = 0;
  T absSoftAppliedFluidPowerW = 0;
  T maxAbsSoftAppliedFluidPowerW = 0;
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
  const std::vector<T>* histMeanSoftFluidTorqueNm = nullptr;
  const std::vector<T>* histMaxSoftFluidTorqueNm = nullptr;
  const std::vector<T>* histMeanSoftAngleStep = nullptr;
  const std::vector<T>* histMaxSoftAngleStep = nullptr;
  const std::vector<T>* histSoftCouplingResidual = nullptr;
  const std::vector<T>* histMaxSoftCouplingResidual = nullptr;
  const std::vector<T>* histSoftCouplingItersUsed = nullptr;
  const std::vector<T>* histMaxSoftCouplingItersUsed = nullptr;
  const std::vector<T>* histSoftElasticEnergyJ = nullptr;
  const std::vector<T>* histMaxSoftElasticEnergyJ = nullptr;
  const std::vector<T>* histSoftDampingPowerW = nullptr;
  const std::vector<T>* histSoftActuatorPowerProxyW = nullptr;
  const std::vector<T>* histAbsSoftActuatorPowerProxyW = nullptr;
  const std::vector<T>* histMaxAbsSoftActuatorPowerProxyW = nullptr;
  const std::vector<T>* histSoftAppliedFluidPowerW = nullptr;
  const std::vector<T>* histAbsSoftAppliedFluidPowerW = nullptr;
  const std::vector<T>* histMaxAbsSoftAppliedFluidPowerW = nullptr;
};

struct SimulationHistory {
  std::vector<T> histT, histVx, histVy, histWz;
  std::vector<T> histFx, histFy, histTz;
  std::vector<T> histXcm, histYcm, histTheta;
  std::vector<T> histPower, histEta;
  std::vector<T> histPowerRigid, histPowerDef;
  std::vector<T> histMeanResidualSlip, histMaxResidualSlip;
  std::vector<T> histMeanSoftFluidTorqueNm, histMaxSoftFluidTorqueNm;
  std::vector<T> histMeanSoftAngleStep, histMaxSoftAngleStep;
  std::vector<T> histSoftCouplingResidual, histMaxSoftCouplingResidual;
  std::vector<T> histSoftCouplingItersUsed, histMaxSoftCouplingItersUsed;
  std::vector<T> histSoftElasticEnergyJ, histMaxSoftElasticEnergyJ;
  std::vector<T> histSoftDampingPowerW;
  std::vector<T> histSoftActuatorPowerProxyW;
  std::vector<T> histAbsSoftActuatorPowerProxyW;
  std::vector<T> histMaxAbsSoftActuatorPowerProxyW;
  std::vector<T> histSoftAppliedFluidPowerW;
  std::vector<T> histAbsSoftAppliedFluidPowerW;
  std::vector<T> histMaxAbsSoftAppliedFluidPowerW;
  std::vector<T> histUswim, histUlat;
  std::vector<T> histForwardNetForce, histFlat;
  std::vector<T> histRe, histSt, histMach, histUstar;
  std::vector<T> histMeanSlip, histMaxSlip;
  std::vector<T> histMeanMarkerForce, histMaxMarkerForce;
  std::vector<T> histNormalizedSlip;

  void append(const HistoryFrameSample& sample);
  HistoryCsvData csvData() const;
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
