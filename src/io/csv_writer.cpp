#include "io/csv_writer.hpp"

#include "io/filesystem.hpp"
#include "physics/material.hpp"
#include "physics/soft_rod.hpp"

#include <fstream>
#include <iomanip>

CsvPathResolution resolveCsvPath(const std::string& requestedPath,
                                 const std::string& expectedHeader,
                                 const std::string& schemaSuffix)
{
  CsvPathResolution result;
  result.path = requestedPath;

  std::ifstream probe(requestedPath);
  if (!probe.good()) {
    return result;
  }

  std::string firstLine;
  if (!std::getline(probe, firstLine) || firstLine.empty()) {
    return result;
  }
  if (firstLine == expectedHeader) {
    return result;
  }

  result.headerMismatch = true;
  result.fallbackPath = appendSuffixBeforeExtension(requestedPath, schemaSuffix);
  result.path = result.fallbackPath;
  return result;
}

CsvWriteResult appendArSummaryCsv(const std::string& requestedPath,
                                   const SummaryCsvInputs& in)
{
  const MaterialProperties material = resolveMaterialProperties(in.p);
  const PlanarRodSectionEstimate rodSection =
    estimatePlanarRodSection(in.p, material);
  const std::string summaryHeader =
    "aspectRatio,bodyAreaTarget,bodyRadius,centerlineLength,totalLength,width,"
    "mass,Ibody,meanU,meanP,meanForwardNetForce,meanLateralForce,"
    "etaNetForceDiagnostic,CoT,meanRe,meanSt,meanUstar,meanEta,"
    "nSteadySamples,simCase,nx,ny,tau,eelFreq,eelLambda,eelA0,"
    "restTime,rampTime,tCut,"
    "meanAbsForwardNetForce,meanSignedForwardNetForce,hydroCost,"
    "transportEfficiencyProxy,nSteadyCycles,cycleMeanUswim,cycleCvUswim,"
    "cycleMeanUstar,cycleCvUstar,"
    "cycleMeanCoT,cycleCvCoT,cycleMeanHydroCost,cycleCvHydroCost,"
    "cycleMeanNormalizedSlip,cycleCvNormalizedSlip,cycleConverged,"
    "warmupMode,waveDirection,initialPlacementClamped,"
    "initialPlacementClampCount,runtimeDomainClampHit,"
    "runtimeDomainClampCount,"
    "alphaIBM,ibmIterations,meanResidualSlip,"
    "maxResidualSlip,meanConstraintPower,meanRigidBodyPower,"
    "meanDeformationPower,meanAbsConstraintPower,meanAbsRigidBodyPower,"
    "meanAbsDeformationPower,transportEfficiencyDef,cotDef,wallBoundary,"
    "bodyKinematics,softBackboneDynamics,"
    "softBackboneRelaxationTime,softBackboneFluidTorqueScale,"
    "softBackboneFluidTorqueFilterTime,"
    "softBackboneAddedMassFrac,softBackboneMaxAngleStep,"
    "softBackboneCouplingIterations,"
    "softBackboneCouplingRelaxation,softBackboneCouplingTolerance,"
    "softBackboneLoadProjection,"
    "meanSoftFluidTorqueNm,"
    "maxSoftFluidTorqueNm,meanSoftAngleStep,maxSoftAngleStep,"
    "meanSoftCouplingResidual,maxSoftCouplingResidual,"
    "meanSoftCouplingItersUsed,maxSoftCouplingItersUsed,"
    "meanSoftElasticEnergyJ,maxSoftElasticEnergyJ,"
    "meanSoftDampingPowerW,meanSoftActuatorPowerProxyW,"
    "meanAbsSoftActuatorPowerProxyW,maxAbsSoftActuatorPowerProxyW,"
    "meanSoftAppliedFluidPowerW,meanAbsSoftAppliedFluidPowerW,"
    "maxAbsSoftAppliedFluidPowerW,softBodyMassKg,meanUPhysicalMps,"
    "cotSoftActuatorProxySI,"
    "bodyMaterial,rhoBodyRatio,youngModulusPa,"
    "poissonRatio,shearModulusPa,bulkModulusPa,materialDampingRatio,"
    "physicalBodyLengthM,"
    "rodWidthM,rodThicknessM,rodBendingStiffnessNm2,rodAxialStiffnessN,"
    "ibmCoupling";

  CsvWriteResult result;
  result.resolution = resolveCsvPath(requestedPath, summaryHeader, "_metrics_v21");

  std::ifstream probe(result.resolution.path);
  const bool writeHeader = !probe.good() || probe.peek() == std::ifstream::traits_type::eof();
  probe.close();

  std::ofstream csv(result.resolution.path, std::ios::app);
  result.opened = csv.is_open();
  if (!result.opened) {
    return result;
  }

  if (writeHeader) {
    csv << summaryHeader << "\n";
  }
  csv << std::setprecision(17)
      << in.p.aspectRatio << ","
      << in.p.bodyAreaTarget << ","
      << in.p.bodyRadius << ","
      << in.bodyCenterlineLength << ","
      << in.bodyLength << ","
      << in.bodyWidth << ","
      << in.mass << ","
      << in.Ibody << ","
      << in.steady.meanU << ","
      << in.steady.meanP << ","
      << in.steady.meanForwardNetForce << ","
      << in.steady.meanLateralForce << ","
      << in.steady.etaNetForceDiagnostic << ","
      << in.steady.CoT << ","
      << in.steady.meanRe << ","
      << in.steady.meanSt << ","
      << in.steady.meanUstar << ","
      << in.steady.meanEta << ","
      << in.steady.nSamples << ","
      << simulationCaseName(in.simCase) << ","
      << in.nx << ","
      << in.ny << ","
      << in.p.tau << ","
      << in.p.eelFreq << ","
      << in.p.eelLambda << ","
      << in.p.eelA0 << ","
      << in.p.restTime << ","
      << in.p.rampTime << ","
      << in.tCut << ","
      << in.steady.meanAbsForwardNetForce << ","
      << in.steady.meanSignedForwardNetForce << ","
      << in.steady.hydroCost << ","
      << in.steady.transportEfficiencyProxy << ","
      << in.cycleConv.nSteadyCycles << ","
      << in.cycleConv.cycleMeanUswim << ","
      << in.cycleConv.cycleCvUswim << ","
      << in.cycleConv.cycleMeanUstar << ","
      << in.cycleConv.cycleCvUstar << ","
      << in.cycleConv.cycleMeanCoT << ","
      << in.cycleConv.cycleCvCoT << ","
      << in.cycleConv.cycleMeanHydroCost << ","
      << in.cycleConv.cycleCvHydroCost << ","
      << in.cycleConv.cycleMeanNormalizedSlip << ","
      << in.cycleConv.cycleCvNormalizedSlip << ","
      << (in.cycleConv.cycleConverged ? 1 : 0) << ","
      << warmupModeName(in.warmupMode) << ","
      << waveDirectionName(in.p.waveDirection) << ","
      << (in.initialPlacementClamped ? 1 : 0) << ","
      << in.initialPlacementClampCount << ","
      << (in.runtimeDomainClampHit ? 1 : 0) << ","
      << in.runtimeDomainClampCount << ","
      << in.alphaIBM << ","
      << in.ibmIterations << ","
      << in.steady.meanResidualSlip << ","
      << in.steady.maxResidualSlip << ","
      << in.steady.meanConstraintPower << ","
      << in.steady.meanRigidBodyPower << ","
      << in.steady.meanDeformationPower << ","
      << in.steady.meanAbsConstraintPower << ","
      << in.steady.meanAbsRigidBodyPower << ","
      << in.steady.meanAbsDeformationPower << ","
      << in.steady.transportEfficiencyDef << ","
      << in.steady.cotDef << ","
      << wallBoundaryName(in.wallBoundary) << ","
      << bodyKinematicsName(in.bodyKinematics) << ","
      << (in.softBackboneDynamics ? 1 : 0) << ","
      << in.softBackboneRelaxationTime << ","
      << in.softBackboneFluidTorqueScale << ","
      << in.softBackboneFluidTorqueFilterTime << ","
      << in.softBackboneAddedMassFrac << ","
      << in.softBackboneMaxAngleStep << ","
      << in.softBackboneCouplingIterations << ","
      << in.softBackboneCouplingRelaxation << ","
      << in.softBackboneCouplingTolerance << ","
      << softBackboneLoadProjectionName(in.softBackboneLoadProjection) << ","
      << in.meanSoftFluidTorqueNm << ","
      << in.maxSoftFluidTorqueNm << ","
      << in.meanSoftAngleStep << ","
      << in.maxSoftAngleStep << ","
      << in.meanSoftCouplingResidual << ","
      << in.maxSoftCouplingResidual << ","
      << in.meanSoftCouplingItersUsed << ","
      << in.maxSoftCouplingItersUsed << ","
      << in.meanSoftElasticEnergyJ << ","
      << in.maxSoftElasticEnergyJ << ","
      << in.meanSoftDampingPowerW << ","
      << in.meanSoftActuatorPowerProxyW << ","
      << in.meanAbsSoftActuatorPowerProxyW << ","
      << in.maxAbsSoftActuatorPowerProxyW << ","
      << in.meanSoftAppliedFluidPowerW << ","
      << in.meanAbsSoftAppliedFluidPowerW << ","
      << in.maxAbsSoftAppliedFluidPowerW << ","
      << in.softBodyMassKg << ","
      << in.meanUPhysicalMps << ","
      << in.cotSoftActuatorProxySI << ","
      << material.name << ","
      << material.densityRatio << ","
      << material.youngModulusPa << ","
      << material.poissonRatio << ","
      << material.shearModulusPa << ","
      << material.bulkModulusPa << ","
      << material.dampingRatio << ","
      << in.p.physicalBodyLengthM << ","
      << (rodSection.valid ? rodSection.widthM : T(0)) << ","
      << (rodSection.valid ? rodSection.thicknessM : T(0)) << ","
      << (rodSection.valid ? rodSection.bendingStiffnessNm2 : T(0)) << ","
      << (rodSection.valid ? rodSection.axialStiffnessN : T(0)) << ","
      << "sparse_eulerian_mdf" << "\n";
  return result;
}

CsvWriteResult appendVerificationSummaryCsv(const std::string& requestedPath,
                                             const VerificationCsvInputs& in)
{
  const MaterialProperties material = resolveMaterialProperties(in.p);
  const PlanarRodSectionEstimate rodSection =
    estimatePlanarRodSection(in.p, material);
  const std::string verificationHeader =
    "run_id,simCase,aspectRatio,bodyAreaTarget,nx,ny,tau,"
    "eelFreq,eelLambda,eelA0,restTime,rampTime,substeps,"
    "dtAnim,dtLbm,Ttotal,spongeWidth,spongeStrength,addedMassFrac,tCut,"
    "meanU,meanP,meanForwardNetForce,meanLateralForce,etaNetForceDiagnostic,"
    "CoT,meanRe,meanSt,meanUstar,meanEta,meanSlip,meanMaxSlip,"
    "meanMarkerForce,meanMaxMarkerForce,meanNormalizedSlip,nSteadySamples,"
    "runtimeDomainClampHit,runtimeDomainClampCount,bodyRadius,centerlineLength,"
    "totalLength,width,mass,Ibody,warnMeanSlip,warnMaxSlip,warnMarkerForce,"
    "studyMode,meanAbsForwardNetForce,meanSignedForwardNetForce,hydroCost,"
    "transportEfficiencyProxy,nSteadyCycles,cycleMeanUswim,cycleCvUswim,"
    "cycleMeanUstar,cycleCvUstar,"
    "cycleMeanCoT,cycleCvCoT,cycleMeanHydroCost,cycleCvHydroCost,"
    "cycleMeanNormalizedSlip,cycleCvNormalizedSlip,cycleConverged,"
    "warmupMode,waveDirection,initialPlacementClamped,"
    "initialPlacementClampCount,"
    "alphaIBM,ibmIterations,meanResidualSlip,"
    "maxResidualSlip,meanConstraintPower,meanRigidBodyPower,"
    "meanDeformationPower,meanAbsConstraintPower,meanAbsRigidBodyPower,"
    "meanAbsDeformationPower,transportEfficiencyDef,cotDef,wallBoundary,"
    "bodyKinematics,softBackboneDynamics,"
    "softBackboneRelaxationTime,softBackboneFluidTorqueScale,"
    "softBackboneFluidTorqueFilterTime,"
    "softBackboneAddedMassFrac,softBackboneMaxAngleStep,"
    "softBackboneCouplingIterations,"
    "softBackboneCouplingRelaxation,softBackboneCouplingTolerance,"
    "softBackboneLoadProjection,"
    "meanSoftFluidTorqueNm,"
    "maxSoftFluidTorqueNm,meanSoftAngleStep,maxSoftAngleStep,"
    "meanSoftCouplingResidual,maxSoftCouplingResidual,"
    "meanSoftCouplingItersUsed,maxSoftCouplingItersUsed,"
    "meanSoftElasticEnergyJ,maxSoftElasticEnergyJ,"
    "meanSoftDampingPowerW,meanSoftActuatorPowerProxyW,"
    "meanAbsSoftActuatorPowerProxyW,maxAbsSoftActuatorPowerProxyW,"
    "meanSoftAppliedFluidPowerW,meanAbsSoftAppliedFluidPowerW,"
    "maxAbsSoftAppliedFluidPowerW,softBodyMassKg,meanUPhysicalMps,"
    "cotSoftActuatorProxySI,"
    "bodyMaterial,rhoBodyRatio,youngModulusPa,"
    "poissonRatio,shearModulusPa,bulkModulusPa,materialDampingRatio,"
    "physicalBodyLengthM,"
    "rodWidthM,rodThicknessM,rodBendingStiffnessNm2,rodAxialStiffnessN,"
    "ibmCoupling";

  CsvWriteResult result;
  result.resolution = resolveCsvPath(requestedPath, verificationHeader, "_metrics_v21");

  std::ifstream probe(result.resolution.path);
  const bool writeHeader = !probe.good() || probe.peek() == std::ifstream::traits_type::eof();
  probe.close();

  std::ofstream csv(result.resolution.path, std::ios::app);
  result.opened = csv.is_open();
  if (!result.opened) {
    return result;
  }

  if (writeHeader) {
    csv << verificationHeader << "\n";
  }
  csv << std::setprecision(17)
      << in.runId << ","
      << simulationCaseName(in.simCase) << ","
      << in.p.aspectRatio << ","
      << in.p.bodyAreaTarget << ","
      << in.nx << ","
      << in.ny << ","
      << in.p.tau << ","
      << in.p.eelFreq << ","
      << in.p.eelLambda << ","
      << in.p.eelA0 << ","
      << in.p.restTime << ","
      << in.p.rampTime << ","
      << in.p.substeps << ","
      << in.p.dtAnim << ","
      << in.dtLbm << ","
      << in.p.Ttotal << ","
      << in.p.spongeWidth << ","
      << in.p.spongeStrength << ","
      << in.p.addedMassFrac << ","
      << in.tCut << ","
      << in.steady.meanU << ","
      << in.steady.meanP << ","
      << in.steady.meanForwardNetForce << ","
      << in.steady.meanLateralForce << ","
      << in.steady.etaNetForceDiagnostic << ","
      << in.steady.CoT << ","
      << in.steady.meanRe << ","
      << in.steady.meanSt << ","
      << in.steady.meanUstar << ","
      << in.steady.meanEta << ","
      << in.steady.meanSlip << ","
      << in.steady.meanMaxSlip << ","
      << in.steady.meanMarkerForce << ","
      << in.steady.meanMaxMarkerForce << ","
      << in.steady.meanNormalizedSlip << ","
      << in.steady.nSamples << ","
      << (in.runtimeDomainClampHit ? 1 : 0) << ","
      << in.runtimeDomainClampCount << ","
      << in.p.bodyRadius << ","
      << in.bodyCenterlineLength << ","
      << in.bodyLength << ","
      << in.bodyWidth << ","
      << in.mass << ","
      << in.Ibody << ","
      << in.p.warnMeanSlip << ","
      << in.p.warnMaxSlip << ","
      << in.p.warnMarkerForce << ","
      << studyModeName(in.studyMode) << ","
      << in.steady.meanAbsForwardNetForce << ","
      << in.steady.meanSignedForwardNetForce << ","
      << in.steady.hydroCost << ","
      << in.steady.transportEfficiencyProxy << ","
      << in.cycleConv.nSteadyCycles << ","
      << in.cycleConv.cycleMeanUswim << ","
      << in.cycleConv.cycleCvUswim << ","
      << in.cycleConv.cycleMeanUstar << ","
      << in.cycleConv.cycleCvUstar << ","
      << in.cycleConv.cycleMeanCoT << ","
      << in.cycleConv.cycleCvCoT << ","
      << in.cycleConv.cycleMeanHydroCost << ","
      << in.cycleConv.cycleCvHydroCost << ","
      << in.cycleConv.cycleMeanNormalizedSlip << ","
      << in.cycleConv.cycleCvNormalizedSlip << ","
      << (in.cycleConv.cycleConverged ? 1 : 0) << ","
      << warmupModeName(in.warmupMode) << ","
      << waveDirectionName(in.p.waveDirection) << ","
      << (in.initialPlacementClamped ? 1 : 0) << ","
      << in.initialPlacementClampCount << ","
      << in.alphaIBM << ","
      << in.ibmIterations << ","
      << in.steady.meanResidualSlip << ","
      << in.steady.maxResidualSlip << ","
      << in.steady.meanConstraintPower << ","
      << in.steady.meanRigidBodyPower << ","
      << in.steady.meanDeformationPower << ","
      << in.steady.meanAbsConstraintPower << ","
      << in.steady.meanAbsRigidBodyPower << ","
      << in.steady.meanAbsDeformationPower << ","
      << in.steady.transportEfficiencyDef << ","
      << in.steady.cotDef << ","
      << wallBoundaryName(in.wallBoundary) << ","
      << bodyKinematicsName(in.bodyKinematics) << ","
      << (in.softBackboneDynamics ? 1 : 0) << ","
      << in.softBackboneRelaxationTime << ","
      << in.softBackboneFluidTorqueScale << ","
      << in.softBackboneFluidTorqueFilterTime << ","
      << in.softBackboneAddedMassFrac << ","
      << in.softBackboneMaxAngleStep << ","
      << in.softBackboneCouplingIterations << ","
      << in.softBackboneCouplingRelaxation << ","
      << in.softBackboneCouplingTolerance << ","
      << softBackboneLoadProjectionName(in.softBackboneLoadProjection) << ","
      << in.meanSoftFluidTorqueNm << ","
      << in.maxSoftFluidTorqueNm << ","
      << in.meanSoftAngleStep << ","
      << in.maxSoftAngleStep << ","
      << in.meanSoftCouplingResidual << ","
      << in.maxSoftCouplingResidual << ","
      << in.meanSoftCouplingItersUsed << ","
      << in.maxSoftCouplingItersUsed << ","
      << in.meanSoftElasticEnergyJ << ","
      << in.maxSoftElasticEnergyJ << ","
      << in.meanSoftDampingPowerW << ","
      << in.meanSoftActuatorPowerProxyW << ","
      << in.meanAbsSoftActuatorPowerProxyW << ","
      << in.maxAbsSoftActuatorPowerProxyW << ","
      << in.meanSoftAppliedFluidPowerW << ","
      << in.meanAbsSoftAppliedFluidPowerW << ","
      << in.maxAbsSoftAppliedFluidPowerW << ","
      << in.softBodyMassKg << ","
      << in.meanUPhysicalMps << ","
      << in.cotSoftActuatorProxySI << ","
      << material.name << ","
      << material.densityRatio << ","
      << material.youngModulusPa << ","
      << material.poissonRatio << ","
      << material.shearModulusPa << ","
      << material.bulkModulusPa << ","
      << material.dampingRatio << ","
      << in.p.physicalBodyLengthM << ","
      << (rodSection.valid ? rodSection.widthM : T(0)) << ","
      << (rodSection.valid ? rodSection.thicknessM : T(0)) << ","
      << (rodSection.valid ? rodSection.bendingStiffnessNm2 : T(0)) << ","
      << (rodSection.valid ? rodSection.axialStiffnessN : T(0)) << ","
      << "sparse_eulerian_mdf" << "\n";
  return result;
}

void SimulationHistory::append(const HistoryFrameSample& s)
{
  histT.push_back(s.t);
  histVx.push_back(s.Vx);
  histVy.push_back(s.Vy);
  histWz.push_back(s.omegaZ);
  histFx.push_back(s.Fx);
  histFy.push_back(s.Fy);
  histTz.push_back(s.Tz);
  histXcm.push_back(s.xCm);
  histYcm.push_back(s.yCm);
  histTheta.push_back(s.theta);
  histPower.push_back(s.power);
  histEta.push_back(s.eta);
  histUswim.push_back(s.Uswim);
  histUlat.push_back(s.Ulateral);
  histForwardNetForce.push_back(s.forwardNetForce);
  histFlat.push_back(s.Flat);
  histRe.push_back(s.Re);
  histSt.push_back(s.St);
  histMach.push_back(s.Mach);
  histUstar.push_back(s.Ustar);
  histMeanSlip.push_back(s.meanSlip);
  histMaxSlip.push_back(s.maxSlip);
  histMeanMarkerForce.push_back(s.meanMarkerForce);
  histMaxMarkerForce.push_back(s.maxMarkerForce);
  histNormalizedSlip.push_back(s.normalizedSlip);
  histPowerRigid.push_back(s.powerRigid);
  histPowerDef.push_back(s.powerDef);
  histMeanResidualSlip.push_back(s.meanResidualSlip);
  histMaxResidualSlip.push_back(s.maxResidualSlip);
  histMeanSoftFluidTorqueNm.push_back(s.meanSoftFluidTorqueNm);
  histMaxSoftFluidTorqueNm.push_back(s.maxSoftFluidTorqueNm);
  histMeanSoftAngleStep.push_back(s.meanSoftAngleStep);
  histMaxSoftAngleStep.push_back(s.maxSoftAngleStep);
  histSoftCouplingResidual.push_back(s.softCouplingResidual);
  histMaxSoftCouplingResidual.push_back(s.maxSoftCouplingResidual);
  histSoftCouplingItersUsed.push_back(s.softCouplingItersUsed);
  histMaxSoftCouplingItersUsed.push_back(s.maxSoftCouplingItersUsed);
  histSoftElasticEnergyJ.push_back(s.softElasticEnergyJ);
  histMaxSoftElasticEnergyJ.push_back(s.maxSoftElasticEnergyJ);
  histSoftDampingPowerW.push_back(s.softDampingPowerW);
  histSoftActuatorPowerProxyW.push_back(s.softActuatorPowerProxyW);
  histAbsSoftActuatorPowerProxyW.push_back(s.absSoftActuatorPowerProxyW);
  histMaxAbsSoftActuatorPowerProxyW.push_back(s.maxAbsSoftActuatorPowerProxyW);
  histSoftAppliedFluidPowerW.push_back(s.softAppliedFluidPowerW);
  histAbsSoftAppliedFluidPowerW.push_back(s.absSoftAppliedFluidPowerW);
  histMaxAbsSoftAppliedFluidPowerW.push_back(s.maxAbsSoftAppliedFluidPowerW);
}

HistoryCsvData SimulationHistory::csvData() const
{
  HistoryCsvData data;
  data.histT = &histT;
  data.histVx = &histVx;
  data.histVy = &histVy;
  data.histWz = &histWz;
  data.histFx = &histFx;
  data.histFy = &histFy;
  data.histTz = &histTz;
  data.histXcm = &histXcm;
  data.histYcm = &histYcm;
  data.histTheta = &histTheta;
  data.histPower = &histPower;
  data.histEta = &histEta;
  data.histUswim = &histUswim;
  data.histUlat = &histUlat;
  data.histForwardNetForce = &histForwardNetForce;
  data.histFlat = &histFlat;
  data.histRe = &histRe;
  data.histSt = &histSt;
  data.histMach = &histMach;
  data.histUstar = &histUstar;
  data.histMeanSlip = &histMeanSlip;
  data.histMaxSlip = &histMaxSlip;
  data.histMeanMarkerForce = &histMeanMarkerForce;
  data.histMaxMarkerForce = &histMaxMarkerForce;
  data.histNormalizedSlip = &histNormalizedSlip;
  data.histPowerRigid = &histPowerRigid;
  data.histPowerDef = &histPowerDef;
  data.histMeanResidualSlip = &histMeanResidualSlip;
  data.histMaxResidualSlip = &histMaxResidualSlip;
  data.histMeanSoftFluidTorqueNm = &histMeanSoftFluidTorqueNm;
  data.histMaxSoftFluidTorqueNm = &histMaxSoftFluidTorqueNm;
  data.histMeanSoftAngleStep = &histMeanSoftAngleStep;
  data.histMaxSoftAngleStep = &histMaxSoftAngleStep;
  data.histSoftCouplingResidual = &histSoftCouplingResidual;
  data.histMaxSoftCouplingResidual = &histMaxSoftCouplingResidual;
  data.histSoftCouplingItersUsed = &histSoftCouplingItersUsed;
  data.histMaxSoftCouplingItersUsed = &histMaxSoftCouplingItersUsed;
  data.histSoftElasticEnergyJ = &histSoftElasticEnergyJ;
  data.histMaxSoftElasticEnergyJ = &histMaxSoftElasticEnergyJ;
  data.histSoftDampingPowerW = &histSoftDampingPowerW;
  data.histSoftActuatorPowerProxyW = &histSoftActuatorPowerProxyW;
  data.histAbsSoftActuatorPowerProxyW = &histAbsSoftActuatorPowerProxyW;
  data.histMaxAbsSoftActuatorPowerProxyW = &histMaxAbsSoftActuatorPowerProxyW;
  data.histSoftAppliedFluidPowerW = &histSoftAppliedFluidPowerW;
  data.histAbsSoftAppliedFluidPowerW = &histAbsSoftAppliedFluidPowerW;
  data.histMaxAbsSoftAppliedFluidPowerW = &histMaxAbsSoftAppliedFluidPowerW;
  return data;
}

std::string writeHistoryCsv(const std::string& runOutDir, const HistoryCsvData& h)
{
  std::string fname = runOutDir + "eel3dof_history.csv";
  std::ofstream csv(fname);
  csv << "t,Vx,Vy,omega_z,Fx,Fy,Tz,x_cm,y_cm,theta,power,eta,"
      << "Uswim,Ulateral,forwardNetForce,Flat,Re,St,Mach,Ustar,"
      << "meanSlip,maxSlip,meanMarkerForce,maxMarkerForce,normalizedSlip,"
      << "constraintPower,rigidBodyPower,deformationPower,"
      << "meanResidualSlip,maxResidualSlip,"
      << "meanSoftFluidTorqueNm,maxSoftFluidTorqueNm,"
      << "meanSoftAngleStep,maxSoftAngleStep,"
      << "softCouplingResidual,maxSoftCouplingResidual,"
      << "softCouplingItersUsed,maxSoftCouplingItersUsed,"
      << "softElasticEnergyJ,maxSoftElasticEnergyJ,"
      << "softDampingPowerW,softActuatorPowerProxyW,"
      << "absSoftActuatorPowerProxyW,maxAbsSoftActuatorPowerProxyW,"
      << "softAppliedFluidPowerW,absSoftAppliedFluidPowerW,"
      << "maxAbsSoftAppliedFluidPowerW\n";
  csv << std::setprecision(17);
  auto atOrZero = [](const std::vector<T>* values, size_t i) -> T {
    return (values != nullptr && i < values->size()) ? (*values)[i] : T(0);
  };
  for (size_t i = 0; i < h.histT->size(); ++i) {
    csv << (*h.histT)[i] << ","
        << (*h.histVx)[i] << "," << (*h.histVy)[i] << "," << (*h.histWz)[i] << ","
        << (*h.histFx)[i] << "," << (*h.histFy)[i] << "," << (*h.histTz)[i] << ","
        << (*h.histXcm)[i] << "," << (*h.histYcm)[i] << "," << (*h.histTheta)[i] << ","
        << (*h.histPower)[i] << "," << (*h.histEta)[i] << ","
        << (*h.histUswim)[i] << "," << (*h.histUlat)[i] << ","
        << (*h.histForwardNetForce)[i] << "," << (*h.histFlat)[i] << ","
        << (*h.histRe)[i] << "," << (*h.histSt)[i] << ","
        << (*h.histMach)[i] << "," << (*h.histUstar)[i] << ","
        << (*h.histMeanSlip)[i] << "," << (*h.histMaxSlip)[i] << ","
        << (*h.histMeanMarkerForce)[i] << "," << (*h.histMaxMarkerForce)[i] << ","
        << (*h.histNormalizedSlip)[i] << ","
        << (*h.histPower)[i] << ","
        << (*h.histPowerRigid)[i] << "," << (*h.histPowerDef)[i] << ","
        << (*h.histMeanResidualSlip)[i] << "," << (*h.histMaxResidualSlip)[i] << ","
        << atOrZero(h.histMeanSoftFluidTorqueNm, i) << ","
        << atOrZero(h.histMaxSoftFluidTorqueNm, i) << ","
        << atOrZero(h.histMeanSoftAngleStep, i) << ","
        << atOrZero(h.histMaxSoftAngleStep, i) << ","
        << atOrZero(h.histSoftCouplingResidual, i) << ","
        << atOrZero(h.histMaxSoftCouplingResidual, i) << ","
        << atOrZero(h.histSoftCouplingItersUsed, i) << ","
        << atOrZero(h.histMaxSoftCouplingItersUsed, i) << ","
        << atOrZero(h.histSoftElasticEnergyJ, i) << ","
        << atOrZero(h.histMaxSoftElasticEnergyJ, i) << ","
        << atOrZero(h.histSoftDampingPowerW, i) << ","
        << atOrZero(h.histSoftActuatorPowerProxyW, i) << ","
        << atOrZero(h.histAbsSoftActuatorPowerProxyW, i) << ","
        << atOrZero(h.histMaxAbsSoftActuatorPowerProxyW, i) << ","
        << atOrZero(h.histSoftAppliedFluidPowerW, i) << ","
        << atOrZero(h.histAbsSoftAppliedFluidPowerW, i) << ","
        << atOrZero(h.histMaxAbsSoftAppliedFluidPowerW, i) << "\n";
  }
  csv.close();
  return fname;
}

std::string writeFinalBodyCsv(const std::string& runOutDir, const LagrangianMarkers& markers)
{
  std::string fname = runOutDir + "eel3dof_body_final.csv";
  std::ofstream csv(fname);
  csv << "x,y\n";
  for (int m = 0; m < markers.size(); ++m) {
    csv << markers.x[m] << "," << markers.y[m] << "\n";
  }
  csv.close();
  return fname;
}
