#include "io/cli.hpp"

#include <cmath>
#include <exception>
#include <sstream>
#include <string>

namespace {

void addError(CliDiagnostics& diagnostics, const std::string& message)
{
  diagnostics.errors.push_back(message);
}

void addWarning(CliDiagnostics& diagnostics, const std::string& message)
{
  diagnostics.warnings.push_back(message);
}

bool parseIntValue(const std::string& key, const std::string& value,
                   int& out, CliDiagnostics& diagnostics)
{
  try {
    size_t pos = 0;
    const int parsed = std::stoi(value, &pos);
    if (pos != value.size()) {
      addError(diagnostics, key + " expects an integer, got '" + value + "'.");
      return false;
    }
    out = parsed;
    return true;
  } catch (const std::exception&) {
    addError(diagnostics, key + " expects an integer, got '" + value + "'.");
    return false;
  }
}

bool parseRealValue(const std::string& key, const std::string& value,
                    T& out, CliDiagnostics& diagnostics)
{
  try {
    size_t pos = 0;
    const T parsed = std::stod(value, &pos);
    if (pos != value.size() || !std::isfinite(parsed)) {
      addError(diagnostics, key + " expects a finite number, got '" + value + "'.");
      return false;
    }
    out = parsed;
    return true;
  } catch (const std::exception&) {
    addError(diagnostics, key + " expects a finite number, got '" + value + "'.");
    return false;
  }
}

bool parseBoolValue(const std::string& key, const std::string& value,
                    bool& out, CliDiagnostics& diagnostics)
{
  if (value == "1" || value == "true" || value == "yes" || value == "on") {
    out = true;
    return true;
  }
  if (value == "0" || value == "false" || value == "no" || value == "off") {
    out = false;
    return true;
  }
  addError(diagnostics, key + " expects a boolean: true/false, 1/0, yes/no, or on/off.");
  return false;
}

bool parseRunModeValue(const std::string& key, const std::string& value,
                       RunMode& out, CliDiagnostics& diagnostics)
{
  if (value == "preview") {
    out = RunMode::Preview;
    return true;
  }
  if (value == "standard") {
    out = RunMode::Standard;
    return true;
  }
  if (value == "full") {
    out = RunMode::Full;
    return true;
  }
  addError(diagnostics, key + " must be preview, standard, or full.");
  return false;
}

bool parseStudyModeValue(const std::string& key, const std::string& value,
                         StudyMode& out, CliDiagnostics& diagnostics)
{
  if (value == "standard") {
    out = StudyMode::Standard;
    return true;
  }
  if (value == "verification" || value == "ar_sweep" || value == "sweep") {
    out = StudyMode::Verification;
    return true;
  }
  addError(diagnostics, key + " must be standard or verification.");
  return false;
}

bool parseCaseValue(const std::string& key, const std::string& value,
                    SimulationCase& out, CliDiagnostics& diagnostics)
{
  if (value == "surge_only") {
    out = SimulationCase::SurgeOnly;
    return true;
  }
  addError(diagnostics,
           key + " only supports surge_only; retired cases are no longer run.");
  return false;
}

bool parseWarmupModeValue(const std::string& key, const std::string& value,
                          WarmupMode& out, CliDiagnostics& diagnostics)
{
  if (value == "rest") {
    out = WarmupMode::Rest;
    return true;
  }
  if (value == "none") {
    out = WarmupMode::None;
    return true;
  }
  addError(diagnostics, key + " must be rest or none.");
  return false;
}

bool parseWallBoundaryValue(const std::string& key, const std::string& value,
                            WallBoundary& out, CliDiagnostics& diagnostics)
{
  if (value == "noslip" || value == "no_slip" || value == "bounceback") {
    out = WallBoundary::NoSlip;
    return true;
  }
  if (value == "freeslip" || value == "free_slip" || value == "fullslip") {
    out = WallBoundary::FreeSlip;
    return true;
  }
  addError(diagnostics, key + " must be noslip or freeslip.");
  return false;
}

bool parseWaveDirectionValue(const std::string& key, const std::string& value,
                             WaveDirection& out, CliDiagnostics& diagnostics)
{
  if (value == "head_to_tail" || value == "headToTail" ||
      value == "head-tail" || value == "forward" || value == "1") {
    out = WaveDirection::HeadToTail;
    return true;
  }
  if (value == "tail_to_head" || value == "tailToHead" ||
      value == "tail-head" || value == "reverse" || value == "-1") {
    out = WaveDirection::TailToHead;
    return true;
  }
  addError(diagnostics, key + " must be head_to_tail or tail_to_head.");
  return false;
}

bool parseBodyKinematicsValue(const std::string& key, const std::string& value,
                              BodyKinematics& out, CliDiagnostics& diagnostics)
{
  if (value == "prescribed_wave" || value == "prescribedWave" ||
      value == "prescribed" || value == "legacy") {
    out = BodyKinematics::PrescribedWave;
    return true;
  }
  if (value == "soft_backbone" || value == "softBackbone" ||
      value == "backbone" || value == "soft") {
    out = BodyKinematics::SoftBackbone;
    return true;
  }
  addError(diagnostics, key + " must be prescribed_wave or soft_backbone.");
  return false;
}

bool parseBodyMaterialValue(const std::string& key, const std::string& value,
                            BodyMaterial& out, CliDiagnostics& diagnostics)
{
  if (value == "neutral" || value == "neutral_lattice" || value == "lattice") {
    out = BodyMaterial::NeutralLattice;
    return true;
  }
  if (value == "dragon_skin_20" || value == "dragonskin20" ||
      value == "dragonSkin20" || value == "ds20") {
    out = BodyMaterial::DragonSkin20;
    return true;
  }
  addError(diagnostics, key + " must be dragon_skin_20 or neutral_lattice.");
  return false;
}

bool parseSoftBackboneLoadProjectionValue(
    const std::string& key, const std::string& value,
    SoftBackboneLoadProjection& out, CliDiagnostics& diagnostics)
{
  if (value == "segment_centroid" || value == "segmentCentroid" ||
      value == "centroid" || value == "legacy" || value == "nearest") {
    out = SoftBackboneLoadProjection::SegmentCentroid;
    return true;
  }
  if (value == "virtual_work" || value == "virtualWork" ||
      value == "cross_section" || value == "crossSection" ||
      value == "cross_section_virtual_work") {
    out = SoftBackboneLoadProjection::CrossSectionVirtualWork;
    return true;
  }
  addError(diagnostics,
           key + " must be segment_centroid or cross_section_virtual_work.");
  return false;
}

void requireFinitePositive(const std::string& name, T value,
                           CliDiagnostics& diagnostics)
{
  if (!std::isfinite(value) || !(value > T(0))) {
    addError(diagnostics, name + " must be > 0.");
  }
}

void requireFiniteNonNegative(const std::string& name, T value,
                              CliDiagnostics& diagnostics)
{
  if (!std::isfinite(value) || value < T(0)) {
    addError(diagnostics, name + " must be >= 0.");
  }
}

}  // namespace

CliParseResult parseCommandLineDetailed(int argc, char* argv[],
                                        const std::string& baseLogDir)
{
  CliParseResult result;
  RunConfig& config = result.config;
  CliDiagnostics& diagnostics = result.diagnostics;
  config.summaryCsv = baseLogDir + "eel3dof_ar_summary.csv";
  config.sensitivityCsv = baseLogDir + "eel3dof_sensitivity_summary.csv";

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "--help" || arg == "-h") {
      result.helpRequested = true;
      continue;
    }

    auto pos = arg.find('=');
    if (pos == std::string::npos) {
      addError(diagnostics, "Expected --key=value, got '" + arg + "'.");
      continue;
    }

    const std::string key = arg.substr(0, pos);
    const std::string val = arg.substr(pos + 1);
    bool boolValue = false;
    int intValue = 0;
    T realValue = T(0);

    if (key == "--nx") {
      if (parseIntValue(key, val, config.p.nx, diagnostics)) { }
    } else if (key == "--ny") {
      if (parseIntValue(key, val, config.p.ny, diagnostics)) { }
    } else if (key == "--tau") {
      if (parseRealValue(key, val, config.p.tau, diagnostics)) { }
    } else if (key == "--Ttotal") {
      if (parseRealValue(key, val, config.p.Ttotal, diagnostics)) { }
    } else if (key == "--dtAnim") {
      if (parseRealValue(key, val, config.p.dtAnim, diagnostics)) { }
    } else if (key == "--substeps") {
      if (parseIntValue(key, val, config.p.substeps, diagnostics)) { }
    } else if (key == "--initialPositionFactor") {
      if (parseRealValue(key, val, config.p.initialPositionFactor, diagnostics)) { }
    } else if (key == "--mode") {
      parseRunModeValue(key, val, config.runMode, diagnostics);
    } else if (key == "--studyMode") {
      parseStudyModeValue(key, val, config.studyMode, diagnostics);
    } else if (key == "--verificationMode") {
      if (parseBoolValue(key, val, boolValue, diagnostics)) {
        config.studyMode = boolValue ? StudyMode::Verification : StudyMode::Standard;
      }
    } else if (key == "--case") {
      parseCaseValue(key, val, config.simCase, diagnostics);
    } else if (key == "--inflowU") {
      if (parseRealValue(key, val, realValue, diagnostics)) {
        addWarning(diagnostics,
                   "--inflowU is ignored: the current surge_only solver uses open pressure boundaries and quiescent initialization.");
      }
    } else if (key == "--tCut") {
      if (parseRealValue(key, val, config.tCut, diagnostics)) { }
    } else if (key == "--summaryCsv") {
      config.summaryCsv = val;
    } else if (key == "--sensitivityCsv") {
      config.sensitivityCsv = val;
    } else if (key == "--runTag") {
      config.runTag = val;
    } else if (key == "--exportVelocity") {
      if (parseBoolValue(key, val, config.exportVelocity, diagnostics)) { }
    } else if (key == "--exportVorticity") {
      if (parseBoolValue(key, val, config.exportVorticity, diagnostics)) { }
    } else if (key == "--exportDiagnostics") {
      if (parseBoolValue(key, val, config.exportDiagnostics, diagnostics)) { }
    } else if (key == "--exportBody") {
      if (parseBoolValue(key, val, config.exportBody, diagnostics)) { }
    } else if (key == "--bodyRadius") {
      if (parseRealValue(key, val, config.p.bodyRadius, diagnostics)) {
        config.rawGeometryOverride = true;
      }
    } else if (key == "--eelScale") {
      if (parseRealValue(key, val, config.p.eelScale, diagnostics)) {
        config.rawGeometryOverride = true;
      }
    } else if (key == "--nSpine") {
      if (parseIntValue(key, val, config.p.nSpine, diagnostics)) { }
    } else if (key == "--eelFreq") {
      if (parseRealValue(key, val, config.p.eelFreq, diagnostics)) { }
    } else if (key == "--eelLambda") {
      if (parseRealValue(key, val, config.p.eelLambda, diagnostics)) { }
    } else if (key == "--eelA0") {
      if (parseRealValue(key, val, config.p.eelA0, diagnostics)) { }
    } else if (key == "--restTime") {
      if (parseRealValue(key, val, config.p.restTime, diagnostics)) { }
    } else if (key == "--rampTime") {
      if (parseRealValue(key, val, config.p.rampTime, diagnostics)) { }
    } else if (key == "--waveDirection" || key == "--gaitDirection") {
      parseWaveDirectionValue(key, val, config.p.waveDirection, diagnostics);
    } else if (key == "--aspectRatio") {
      if (parseRealValue(key, val, config.p.aspectRatio, diagnostics)) {
        config.aspectGeometryOverride = true;
      }
    } else if (key == "--bodyAreaTarget") {
      if (parseRealValue(key, val, config.p.bodyAreaTarget, diagnostics)) {
        config.aspectGeometryOverride = true;
      }
    } else if (key == "--useAspectRatioGeometry") {
      if (parseBoolValue(key, val, config.p.useAspectRatioGeometry, diagnostics)) {
        config.aspectGeometryOverride = true;
      }
    } else if (key == "--alphaIBM") {
      if (parseRealValue(key, val, config.alphaIBM, diagnostics)) { }
    } else if (key == "--ibmIterations") {
      if (parseIntValue(key, val, config.ibmIterations, diagnostics)) { }
    } else if (key == "--nIbmIters") {
      if (parseIntValue(key, val, config.legacyNIbmIters, diagnostics)) {
        addWarning(diagnostics,
                   "--nIbmIters is retained only for legacy input compatibility; use --ibmIterations for sparse Eulerian MDF iterations.");
      }
    } else if (key == "--warnMeanSlip") {
      if (parseRealValue(key, val, config.p.warnMeanSlip, diagnostics)) { }
    } else if (key == "--warnMaxSlip") {
      if (parseRealValue(key, val, config.p.warnMaxSlip, diagnostics)) { }
    } else if (key == "--warnMarkerForce") {
      if (parseRealValue(key, val, config.p.warnMarkerForce, diagnostics)) { }
    } else if (key == "--addedMassFrac") {
      if (parseRealValue(key, val, config.p.addedMassFrac, diagnostics)) { }
    } else if (key == "--material" || key == "--bodyMaterial") {
      parseBodyMaterialValue(key, val, config.p.bodyMaterial, diagnostics);
    } else if (key == "--fluidDensityKgM3") {
      if (parseRealValue(key, val, config.p.fluidDensityKgM3, diagnostics)) { }
    } else if (key == "--rhoBodyRatio") {
      if (parseRealValue(key, val, config.p.rhoBodyRatioOverride, diagnostics)) { }
    } else if (key == "--youngModulusPa") {
      if (parseRealValue(key, val, config.p.youngModulusPaOverride, diagnostics)) { }
    } else if (key == "--poissonRatio") {
      if (parseRealValue(key, val, config.p.poissonRatioOverride, diagnostics)) { }
    } else if (key == "--materialDampingRatio") {
      if (parseRealValue(key, val, config.p.materialDampingRatio, diagnostics)) { }
    } else if (key == "--physicalBodyLengthM") {
      if (parseRealValue(key, val, config.p.physicalBodyLengthM, diagnostics)) { }
    } else if (key == "--bodyThicknessM") {
      if (parseRealValue(key, val, config.p.bodyThicknessM, diagnostics)) { }
    } else if (key == "--spongeWidth") {
      if (parseIntValue(key, val, config.p.spongeWidth, diagnostics)) { }
    } else if (key == "--spongeStrength") {
      if (parseRealValue(key, val, config.p.spongeStrength, diagnostics)) { }
    } else if (key == "--nWarmup") {
      if (parseIntValue(key, val, config.p.nWarmup, diagnostics)) { }
    } else if (key == "--warmupMode") {
      parseWarmupModeValue(key, val, config.warmupMode, diagnostics);
    } else if (key == "--bodyKinematics" || key == "--deformationSource") {
      parseBodyKinematicsValue(key, val, config.bodyKinematics, diagnostics);
    } else if (key == "--useSoftBackbone") {
      if (parseBoolValue(key, val, boolValue, diagnostics)) {
        config.bodyKinematics = boolValue
                              ? BodyKinematics::SoftBackbone
                              : BodyKinematics::PrescribedWave;
      }
    } else if (key == "--softBackboneDynamics") {
      if (parseBoolValue(key, val, config.softBackboneDynamics, diagnostics)) { }
    } else if (key == "--softBackboneRelaxationTime") {
      if (parseRealValue(key, val, config.softBackboneRelaxationTime, diagnostics)) {
        addWarning(diagnostics,
                   "--softBackboneRelaxationTime is currently reported for schema compatibility but is not used by the implicit backbone integrator.");
      }
    } else if (key == "--softBackboneFluidTorqueScale") {
      if (parseRealValue(key, val, config.softBackboneFluidTorqueScale, diagnostics)) { }
    } else if (key == "--softBackboneFluidTorqueFilterTime") {
      if (parseRealValue(key, val, config.softBackboneFluidTorqueFilterTime, diagnostics)) { }
    } else if (key == "--softBackboneMaxAngleStep") {
      if (parseRealValue(key, val, config.softBackboneMaxAngleStep, diagnostics)) { }
    } else if (key == "--softBackboneCouplingIterations") {
      if (parseIntValue(key, val, config.softBackboneCouplingIterations, diagnostics)) { }
    } else if (key == "--softBackboneCouplingRelaxation") {
      if (parseRealValue(key, val, config.softBackboneCouplingRelaxation, diagnostics)) { }
    } else if (key == "--softBackboneCouplingTolerance") {
      if (parseRealValue(key, val, config.softBackboneCouplingTolerance, diagnostics)) { }
    } else if (key == "--softBackboneLoadProjection") {
      parseSoftBackboneLoadProjectionValue(
        key, val, config.softBackboneLoadProjection, diagnostics);
    } else if (key == "--softBackboneAddedMassFrac") {
      if (parseRealValue(key, val, config.softBackboneAddedMassFrac, diagnostics)) { }
    } else if (key == "--softBackboneAbortOnInstability") {
      if (parseBoolValue(key, val, config.softBackboneAbortOnInstability, diagnostics)) { }
    } else if (key == "--softBackboneAbortMeanSlip") {
      if (parseRealValue(key, val, config.softBackboneAbortMeanSlip, diagnostics)) { }
    } else if (key == "--softBackboneAbortMaxSlip") {
      if (parseRealValue(key, val, config.softBackboneAbortMaxSlip, diagnostics)) { }
    } else if (key == "--softBackboneAbortSaturatedFrames") {
      if (parseIntValue(key, val, intValue, diagnostics)) {
        config.softBackboneAbortSaturatedFrames = intValue;
      }
    } else if (key == "--wallBoundary") {
      parseWallBoundaryValue(key, val, config.wallBoundary, diagnostics);
    } else {
      addError(diagnostics, "Unknown option '" + key + "'. Use --help for supported options.");
    }
  }

  return result;
}

RunConfig parseCommandLine(int argc, char* argv[], const std::string& baseLogDir)
{
  return parseCommandLineDetailed(argc, argv, baseLogDir).config;
}

CliDiagnostics validateRunConfig(const RunConfig& config)
{
  CliDiagnostics diagnostics;
  const EelParams& p = config.p;

  if (p.nx < 8) addError(diagnostics, "--nx must be at least 8.");
  if (p.ny < 8) addError(diagnostics, "--ny must be at least 8.");
  if (!std::isfinite(p.tau) || !(p.tau > T(0.5))) {
    addError(diagnostics, "--tau must be > 0.5 so the lattice viscosity is positive.");
  }
  requireFinitePositive("--dtAnim", p.dtAnim, diagnostics);
  if (p.substeps < 1) addError(diagnostics, "--substeps must be at least 1.");
  requireFiniteNonNegative("--Ttotal", p.Ttotal, diagnostics);
  if (p.nSpine < 2) addError(diagnostics, "--nSpine must be at least 2.");
  requireFinitePositive("--bodyRadius", p.bodyRadius, diagnostics);
  requireFinitePositive("--eelLength", p.eelLength, diagnostics);
  requireFinitePositive("--eelScale", p.eelScale, diagnostics);
  requireFinitePositive("--eelLambda", p.eelLambda, diagnostics);
  requireFiniteNonNegative("--eelFreq", p.eelFreq, diagnostics);
  requireFiniteNonNegative("--eelA0", p.eelA0, diagnostics);
  requireFiniteNonNegative("--restTime", p.restTime, diagnostics);
  requireFiniteNonNegative("--rampTime", p.rampTime, diagnostics);
  requireFiniteNonNegative("--initialPositionFactor", p.initialPositionFactor, diagnostics);
  requireFiniteNonNegative("--addedMassFrac", p.addedMassFrac, diagnostics);
  requireFiniteNonNegative("--warnMeanSlip", p.warnMeanSlip, diagnostics);
  requireFiniteNonNegative("--warnMaxSlip", p.warnMaxSlip, diagnostics);
  requireFiniteNonNegative("--warnMarkerForce", p.warnMarkerForce, diagnostics);
  requireFinitePositive("--fluidDensityKgM3", p.fluidDensityKgM3, diagnostics);
  if (p.rhoBodyRatioOverride != T(-1) && p.rhoBodyRatioOverride <= T(0)) {
    addError(diagnostics, "--rhoBodyRatio must be > 0 when provided.");
  }
  if (p.youngModulusPaOverride != T(-1) && p.youngModulusPaOverride <= T(0)) {
    addError(diagnostics, "--youngModulusPa must be > 0 when provided.");
  }
  if (p.poissonRatioOverride != T(-1) &&
      !(p.poissonRatioOverride > T(-0.99) && p.poissonRatioOverride < T(0.5))) {
    addError(diagnostics, "--poissonRatio must be in (-0.99, 0.5) when provided.");
  }
  requireFiniteNonNegative("--materialDampingRatio", p.materialDampingRatio, diagnostics);
  requireFiniteNonNegative("--physicalBodyLengthM", p.physicalBodyLengthM, diagnostics);
  requireFiniteNonNegative("--bodyThicknessM", p.bodyThicknessM, diagnostics);
  if (p.spongeWidth < 0) addError(diagnostics, "--spongeWidth must be >= 0.");
  requireFiniteNonNegative("--spongeStrength", p.spongeStrength, diagnostics);
  if (p.nWarmup < 0) addError(diagnostics, "--nWarmup must be >= 0.");
  if (config.tCut < T(0) && config.tCut != T(-1)) {
    addError(diagnostics, "--tCut must be >= 0, or -1 to request the automatic steady-window cut.");
  }
  if (p.useAspectRatioGeometry) {
    if (!std::isfinite(p.aspectRatio) || !(p.aspectRatio > T(1))) {
      addError(diagnostics, "--aspectRatio must be > 1 when aspect-ratio geometry is enabled.");
    }
    requireFinitePositive("--bodyAreaTarget", p.bodyAreaTarget, diagnostics);
  }
  requireFiniteNonNegative("--alphaIBM", config.alphaIBM, diagnostics);
  if (config.ibmIterations < 1) {
    addError(diagnostics, "--ibmIterations must be at least 1.");
  }
  if (config.legacyNIbmIters < 1) {
    addError(diagnostics, "--nIbmIters must be at least 1 when provided.");
  }
  if (config.softBackboneDynamics &&
      config.bodyKinematics != BodyKinematics::SoftBackbone) {
    addError(diagnostics,
             "--softBackboneDynamics=true requires --bodyKinematics=soft_backbone.");
  }
  requireFiniteNonNegative("--softBackboneRelaxationTime",
                           config.softBackboneRelaxationTime, diagnostics);
  requireFiniteNonNegative("--softBackboneFluidTorqueScale",
                           config.softBackboneFluidTorqueScale, diagnostics);
  requireFiniteNonNegative("--softBackboneFluidTorqueFilterTime",
                           config.softBackboneFluidTorqueFilterTime, diagnostics);
  requireFinitePositive("--softBackboneMaxAngleStep",
                        config.softBackboneMaxAngleStep, diagnostics);
  if (config.softBackboneCouplingIterations < 1) {
    addError(diagnostics, "--softBackboneCouplingIterations must be at least 1.");
  }
  if (!std::isfinite(config.softBackboneCouplingRelaxation) ||
      !(config.softBackboneCouplingRelaxation > T(0)) ||
      config.softBackboneCouplingRelaxation > T(1)) {
    addError(diagnostics, "--softBackboneCouplingRelaxation must be in (0, 1].");
  }
  requireFiniteNonNegative("--softBackboneCouplingTolerance",
                           config.softBackboneCouplingTolerance, diagnostics);
  if (!config.softBackboneDynamics &&
      config.softBackboneCouplingIterations > 1) {
    addError(diagnostics,
             "--softBackboneCouplingIterations > 1 requires --softBackboneDynamics=true.");
  }
  requireFiniteNonNegative("--softBackboneAddedMassFrac",
                           config.softBackboneAddedMassFrac, diagnostics);
  requireFiniteNonNegative("--softBackboneAbortMeanSlip",
                           config.softBackboneAbortMeanSlip, diagnostics);
  requireFiniteNonNegative("--softBackboneAbortMaxSlip",
                           config.softBackboneAbortMaxSlip, diagnostics);
  if (config.softBackboneAbortSaturatedFrames < 1) {
    addError(diagnostics, "--softBackboneAbortSaturatedFrames must be at least 1.");
  }
  if (config.summaryCsv.empty()) {
    addError(diagnostics, "--summaryCsv must not be empty.");
  }
  if (config.sensitivityCsv.empty()) {
    addError(diagnostics, "--sensitivityCsv must not be empty.");
  }

  return diagnostics;
}

std::string commandLineUsage(const char* executableName)
{
  std::ostringstream os;
  os << "Usage: " << (executableName ? executableName : "11_lbm_eel_3dof")
     << " [--key=value ...]\n"
     << "\nCommon options:\n"
     << "  --nx=<int> --ny=<int> --tau=<number>\n"
     << "  --Ttotal=<seconds> --dtAnim=<seconds> --substeps=<int>\n"
     << "  --mode=preview|standard|full --studyMode=standard|verification\n"
     << "  --useAspectRatioGeometry=true|false --aspectRatio=<number> --bodyAreaTarget=<number>\n"
     << "  --bodyRadius=<lu> --eelScale=<lu> --nSpine=<int>\n"
     << "  --wallBoundary=noslip|freeslip --waveDirection=head_to_tail|tail_to_head\n"
     << "  --bodyKinematics=prescribed_wave|soft_backbone --softBackboneDynamics=true|false\n"
     << "  --softBackboneFluidTorqueFilterTime=<seconds> --softBackboneCouplingIterations=<int> --softBackboneCouplingTolerance=<rad> --softBackboneLoadProjection=segment_centroid|cross_section_virtual_work\n"
     << "  --exportVelocity=true|false --exportVorticity=true|false --exportDiagnostics=true|false --exportBody=true|false\n"
     << "  --ibmIterations=<int> --alphaIBM=<number>\n"
     << "  --summaryCsv=<path> --sensitivityCsv=<path> --runTag=<token>\n"
     << "\nCompatibility notes:\n"
     << "  --case only supports surge_only.\n"
     << "  --inflowU is ignored by the current open-pressure-boundary solver.\n"
     << "  --nIbmIters is legacy; use --ibmIterations for MDF iterations.\n"
     << "  --softBackboneRelaxationTime is kept only for old CSV schemas.\n"
     << "  Spatial exports default to false; use scripts/run_visualization_case.py for ParaView output.\n";
  return os.str();
}
