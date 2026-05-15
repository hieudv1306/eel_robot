#include "io/cli.hpp"

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

namespace {

CliParseResult parseArgs(const std::vector<std::string>& args)
{
  std::vector<char*> argv;
  argv.reserve(args.size());
  for (const std::string& arg : args) {
    argv.push_back(const_cast<char*>(arg.c_str()));
  }
  return parseCommandLineDetailed(static_cast<int>(argv.size()), argv.data(),
                                  "tmp/");
}

bool containsText(const std::vector<std::string>& values,
                  const std::string& needle)
{
  return std::any_of(values.begin(), values.end(),
                     [&](const std::string& value) {
                       return value.find(needle) != std::string::npos;
                     });
}

}  // namespace

int main()
{
  {
    auto parsed = parseArgs({
      "eel",
      "--nx=320",
      "--ny=120",
      "--tau=0.55",
      "--Ttotal=0.08",
      "--dtAnim=0.04",
      "--substeps=2",
      "--nWarmup=0",
      "--case=surge_only",
      "--mode=preview",
      "--studyMode=verification",
      "--wallBoundary=freeslip",
      "--bodyKinematics=soft_backbone",
      "--softBackboneDynamics=true",
      "--softBackboneCouplingIterations=3",
      "--softBackboneCouplingRelaxation=0.5",
      "--softBackboneCouplingTolerance=1e-5",
      "--softBackboneFluidTorqueFilterTime=0.02",
      "--softBackboneLoadProjection=cross_section_virtual_work"
    });
    assert(parsed.diagnostics.errors.empty());
    const CliDiagnostics validation = validateRunConfig(parsed.config);
    assert(validation.errors.empty());
    assert(parsed.config.p.nx == 320);
    assert(parsed.config.p.ny == 120);
    assert(parsed.config.runMode == RunMode::Preview);
    assert(parsed.config.studyMode == StudyMode::Verification);
    assert(parsed.config.wallBoundary == WallBoundary::FreeSlip);
    assert(!parsed.config.exportVelocity);
    assert(!parsed.config.exportVorticity);
    assert(!parsed.config.exportDiagnostics);
    assert(!parsed.config.exportBody);
    assert(parsed.config.bodyKinematics == BodyKinematics::SoftBackbone);
    assert(parsed.config.p.waveDirection == WaveDirection::HeadToTail);
    assert(parsed.config.softBackboneDynamics);
    assert(parsed.config.softBackboneCouplingIterations == 3);
    assert(parsed.config.softBackboneCouplingRelaxation == 0.5);
    assert(parsed.config.softBackboneCouplingTolerance == 1e-5);
    assert(parsed.config.softBackboneFluidTorqueFilterTime == 0.02);
    assert(parsed.config.softBackboneLoadProjection ==
           SoftBackboneLoadProjection::CrossSectionVirtualWork);
  }

  {
    auto parsed = parseArgs({
      "eel",
      "--exportVelocity=true",
      "--exportVorticity=true",
      "--exportDiagnostics=true",
      "--exportBody=true"
    });
    assert(parsed.diagnostics.errors.empty());
    assert(parsed.config.exportVelocity);
    assert(parsed.config.exportVorticity);
    assert(parsed.config.exportDiagnostics);
    assert(parsed.config.exportBody);
  }

  {
    auto parsed = parseArgs({
      "eel",
      "--bodyKinematics=soft_backbone",
      "--waveDirection=tail_to_head"
    });
    assert(parsed.diagnostics.errors.empty());
    assert(parsed.config.p.waveDirection == WaveDirection::TailToHead);
  }

  {
    auto parsed = parseArgs({"eel", "--unknown=1"});
    assert(!parsed.diagnostics.errors.empty());
    assert(containsText(parsed.diagnostics.errors, "Unknown option"));
  }

  {
    auto parsed = parseArgs({"eel", "--tau=abc"});
    assert(!parsed.diagnostics.errors.empty());
    assert(containsText(parsed.diagnostics.errors, "finite number"));
  }

  {
    auto parsed = parseArgs({"eel", "--case=full3dof"});
    assert(!parsed.diagnostics.errors.empty());
    assert(containsText(parsed.diagnostics.errors, "only supports surge_only"));
  }

  {
    auto parsed = parseArgs({"eel"});
    assert(parsed.config.bodyKinematics == BodyKinematics::SoftBackbone);
  }

  {
    auto parsed = parseArgs({
      "eel",
      "--bodyKinematics=prescribed_wave",
      "--softBackboneDynamics=true"
    });
    const CliDiagnostics validation = validateRunConfig(parsed.config);
    assert(!validation.errors.empty());
    assert(containsText(validation.errors, "requires --bodyKinematics=soft_backbone"));
  }

  {
    auto parsed = parseArgs({
      "eel",
      "--bodyKinematics=soft_backbone",
      "--softBackboneCouplingIterations=2"
    });
    const CliDiagnostics validation = validateRunConfig(parsed.config);
    assert(!validation.errors.empty());
    assert(containsText(validation.errors,
                        "requires --softBackboneDynamics=true"));
  }

  {
    auto parsed = parseArgs({"eel", "--softBackboneCouplingTolerance=-1"});
    const CliDiagnostics validation = validateRunConfig(parsed.config);
    assert(!validation.errors.empty());
    assert(containsText(validation.errors,
                        "--softBackboneCouplingTolerance must be >= 0"));
  }

  {
    auto parsed = parseArgs({"eel", "--softBackboneFluidTorqueFilterTime=-0.1"});
    const CliDiagnostics validation = validateRunConfig(parsed.config);
    assert(!validation.errors.empty());
    assert(containsText(validation.errors,
                        "--softBackboneFluidTorqueFilterTime must be >= 0"));
  }

  {
    auto parsed = parseArgs({
      "eel",
      "--inflowU=0.1",
      "--nIbmIters=2",
      "--softBackboneRelaxationTime=0.2"
    });
    assert(parsed.diagnostics.errors.empty());
    assert(containsText(parsed.diagnostics.warnings, "--inflowU is ignored"));
    assert(containsText(parsed.diagnostics.warnings, "--nIbmIters is retained"));
    assert(containsText(parsed.diagnostics.warnings,
                        "--softBackboneRelaxationTime is currently reported"));
    assert(parsed.config.legacyNIbmIters == 2);
  }

  {
    auto parsed = parseArgs({"eel", "--help"});
    assert(parsed.helpRequested);
    assert(parsed.diagnostics.errors.empty());
    assert(commandLineUsage("eel").find("--nx=<int>") != std::string::npos);
  }
}
