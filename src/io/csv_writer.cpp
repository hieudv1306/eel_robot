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
    "nSteadySamples,simCase,nx,ny,tau,eelFreq,eelLambda,eelA0,tCut,"
    "meanAbsForwardNetForce,meanSignedForwardNetForce,hydroCost,"
    "transportEfficiencyProxy,nSteadyCycles,cycleMeanUstar,cycleCvUstar,"
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
    "softBackboneMaxAngleStep,meanSoftFluidTorqueNm,"
    "maxSoftFluidTorqueNm,meanSoftAngleStep,maxSoftAngleStep,"
    "bodyMaterial,rhoBodyRatio,youngModulusPa,"
    "poissonRatio,shearModulusPa,bulkModulusPa,physicalBodyLengthM,"
    "rodWidthM,rodThicknessM,rodBendingStiffnessNm2,rodAxialStiffnessN,"
    "ibmCoupling";

  CsvWriteResult result;
  result.resolution = resolveCsvPath(requestedPath, summaryHeader, "_metrics_v14");

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
      << in.tCut << ","
      << in.steady.meanAbsForwardNetForce << ","
      << in.steady.meanSignedForwardNetForce << ","
      << in.steady.hydroCost << ","
      << in.steady.transportEfficiencyProxy << ","
      << in.cycleConv.nSteadyCycles << ","
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
      << in.softBackboneMaxAngleStep << ","
      << in.meanSoftFluidTorqueNm << ","
      << in.maxSoftFluidTorqueNm << ","
      << in.meanSoftAngleStep << ","
      << in.maxSoftAngleStep << ","
      << material.name << ","
      << material.densityRatio << ","
      << material.youngModulusPa << ","
      << material.poissonRatio << ","
      << material.shearModulusPa << ","
      << material.bulkModulusPa << ","
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
    "run_id,simCase,aspectRatio,bodyAreaTarget,nx,ny,tau,substeps,"
    "dtAnim,dtLbm,Ttotal,spongeWidth,spongeStrength,addedMassFrac,tCut,"
    "meanU,meanP,meanForwardNetForce,meanLateralForce,etaNetForceDiagnostic,"
    "CoT,meanRe,meanSt,meanUstar,meanEta,meanSlip,meanMaxSlip,"
    "meanMarkerForce,meanMaxMarkerForce,meanNormalizedSlip,nSteadySamples,"
    "runtimeDomainClampHit,runtimeDomainClampCount,bodyRadius,centerlineLength,"
    "totalLength,width,mass,Ibody,warnMeanSlip,warnMaxSlip,warnMarkerForce,"
    "studyMode,meanAbsForwardNetForce,meanSignedForwardNetForce,hydroCost,"
    "transportEfficiencyProxy,nSteadyCycles,cycleMeanUstar,cycleCvUstar,"
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
    "softBackboneMaxAngleStep,meanSoftFluidTorqueNm,"
    "maxSoftFluidTorqueNm,meanSoftAngleStep,maxSoftAngleStep,"
    "bodyMaterial,rhoBodyRatio,youngModulusPa,"
    "poissonRatio,shearModulusPa,bulkModulusPa,physicalBodyLengthM,"
    "rodWidthM,rodThicknessM,rodBendingStiffnessNm2,rodAxialStiffnessN,"
    "ibmCoupling";

  CsvWriteResult result;
  result.resolution = resolveCsvPath(requestedPath, verificationHeader, "_metrics_v14");

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
      << in.softBackboneMaxAngleStep << ","
      << in.meanSoftFluidTorqueNm << ","
      << in.maxSoftFluidTorqueNm << ","
      << in.meanSoftAngleStep << ","
      << in.maxSoftAngleStep << ","
      << material.name << ","
      << material.densityRatio << ","
      << material.youngModulusPa << ","
      << material.poissonRatio << ","
      << material.shearModulusPa << ","
      << material.bulkModulusPa << ","
      << in.p.physicalBodyLengthM << ","
      << (rodSection.valid ? rodSection.widthM : T(0)) << ","
      << (rodSection.valid ? rodSection.thicknessM : T(0)) << ","
      << (rodSection.valid ? rodSection.bendingStiffnessNm2 : T(0)) << ","
      << (rodSection.valid ? rodSection.axialStiffnessN : T(0)) << ","
      << "sparse_eulerian_mdf" << "\n";
  return result;
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
      << "meanSoftAngleStep,maxSoftAngleStep\n";
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
        << atOrZero(h.histMaxSoftAngleStep, i) << "\n";
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
