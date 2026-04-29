#include "io/cli.hpp"

#include <algorithm>
#include <string>

RunConfig parseCommandLine(int argc, char* argv[], const std::string& baseLogDir)
{
  RunConfig config;
  config.summaryCsv = baseLogDir + "eel3dof_ar_summary.csv";
  config.sensitivityCsv = baseLogDir + "eel3dof_kappa_sensitivity_summary.csv";

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    auto pos = arg.find('=');
    if (pos == std::string::npos) continue;
    std::string key = arg.substr(0, pos);
    std::string val = arg.substr(pos + 1);
    if (key == "--nx")            config.p.nx = std::stoi(val);
    if (key == "--ny")            config.p.ny = std::stoi(val);
    if (key == "--tau")           config.p.tau = std::stod(val);
    if (key == "--Ttotal")        config.p.Ttotal = std::stod(val);
    if (key == "--substeps")      config.p.substeps = std::stoi(val);
    if (key == "--initialPositionFactor") config.p.initialPositionFactor = std::stod(val);
    if (key == "--mode")          config.runMode = parseRunMode(val);
    if (key == "--studyMode")     config.studyMode = parseStudyMode(val);
    if (key == "--verificationMode") {
      config.studyMode = parseBool(val) ? StudyMode::Verification : StudyMode::Standard;
    }
    if (key == "--case")          config.simCase = parseSimulationCase(val);
    if (key == "--inflowU")       config.p.inflowVelocity = std::stod(val);
    if (key == "--tCut")          config.tCut = std::stod(val);
    if (key == "--summaryCsv")    config.summaryCsv = val;
    if (key == "--sensitivityCsv") config.sensitivityCsv = val;
    if (key == "--runTag")        config.runTag = val;
    if (key == "--bodyRadius")   { config.p.bodyRadius = std::stod(val); config.rawGeometryOverride = true; }
    if (key == "--eelScale")     { config.p.eelScale = std::stod(val); config.rawGeometryOverride = true; }
    if (key == "--nSpine")        config.p.nSpine = std::stoi(val);
    if (key == "--eelFreq")       config.p.eelFreq = std::stod(val);
    if (key == "--eelLambda")     config.p.eelLambda = std::stod(val);
    if (key == "--eelA0")         config.p.eelA0 = std::stod(val);
    if (key == "--aspectRatio")  { config.p.aspectRatio = std::stod(val); config.aspectGeometryOverride = true; }
    if (key == "--bodyAreaTarget") {
      config.p.bodyAreaTarget = std::stod(val);
      config.aspectGeometryOverride = true;
    }
    if (key == "--useAspectRatioGeometry") {
      config.p.useAspectRatioGeometry = parseBool(val);
      config.aspectGeometryOverride = true;
    }
    if (key == "--kappa") {
      config.p.kappa = std::stod(val);
      config.legacyKappaInputUsed = true;
      config.legacyKappaInputValue = config.p.kappa;
      config.alphaIBM = config.p.kappa;
    }
    if (key == "--alphaIBM")      config.alphaIBM = std::stod(val);
    if (key == "--ibmIterations") config.ibmIterations = std::max(1, std::stoi(val));
    if (key == "--nIbmIters")     config.p.nIbmIters = std::stoi(val);
    if (key == "--warnMeanSlip")  config.p.warnMeanSlip = std::stod(val);
    if (key == "--warnMaxSlip")   config.p.warnMaxSlip = std::stod(val);
    if (key == "--warnMarkerForce") config.p.warnMarkerForce = std::stod(val);
    if (key == "--addedMassFrac") config.p.addedMassFrac = std::stod(val);
    if (key == "--spongeWidth")   config.p.spongeWidth = std::stoi(val);
    if (key == "--spongeStrength")config.p.spongeStrength = std::stod(val);
    if (key == "--nWarmup")       config.p.nWarmup = std::stoi(val);
    if (key == "--warmupMode")    config.warmupMode = parseWarmupMode(val);
    if (key == "--gaitNormalization") config.gaitNormalization = parseGaitNormalization(val);
    if (key == "--tailAmpRatioTarget") config.tailAmpRatioTarget = std::stod(val);
    if (key == "--targetSt")      config.targetSt = std::stod(val);
    if (key == "--referenceU")    config.referenceU = std::stod(val);
  }

  return config;
}
