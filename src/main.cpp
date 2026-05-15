/*
 * Self-propelled eel robot — 3-DOF free swimming (surge, heave, yaw).
 * LBM-IBM with Peskin 4-point delta, Guo forcing, direct-forcing IBM.
 * Newton-Euler equations update (Vx, Vy, omega_z) from IBM forces/torques.
 *
 * Ported from 11_lbm_eel_3dof.py (Taichi) to OpenLB (C++).
 *
 * ================================================================
 *  UNIT CONVENTION (lattice Boltzmann natural units)
 * ================================================================
 *  All quantities inside the per-step loop use LBM natural units:
 *    δx = 1 lu   (lattice spacing)
 *    δt = 1 step (one LBM collide-stream cycle)
 *
 *  Velocities : lu / step   (Vx, Vy, ux, uy, marker Ud/Vd)
 *  Angular vel: rad / step   (omegaZ)
 *  Force      : lu² / step   (IBM body-reaction force, for 2-D)
 *  Mass       : lu²           (area × density, for 2-D neutrally buoyant)
 *  Accel      : lu / step²   (F / M)
 *
 *  Newton-Euler integration per step (δt = 1):
 *    ΔV = F / M              (no extra dt factor)
 *    Δx = V                  (no extra dt factor)
 *
 *  Body-wave kinematics (frequency, ramp rate) are defined in an
 *  external "physical-time" coordinate t [seconds].  The mapping
 *  factor  dtLbm = dtAnim / substeps  [seconds / step]  converts
 *  the analytical ∂y/∂t (lu / s) to the per-step deformation
 *  velocity (lu / step) used by the IBM markers.
 * ================================================================
 */

#include <olb.h>
#include <cmath>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <numeric>
#include <sstream>
#include <cctype>
#include <omp.h>

#include "core/enums.hpp"
#include "core/params.hpp"
#include "core/types.hpp"
#include "io/cli.hpp"
#include "io/filesystem.hpp"
#include "io/csv_writer.hpp"
#include "physics/geometry.hpp"
#include "physics/markers.hpp"
#include "physics/material.hpp"
#include "physics/diagnostics.hpp"
#include "physics/rigid_body.hpp"
#include "physics/soft_backbone.hpp"
#include "physics/soft_rod.hpp"
#include "solver/boundary.hpp"
#include "solver/export_vtk.hpp"
#include "solver/ibm.hpp"

using namespace olb;
using namespace olb::descriptors;


// ============================================================
//  Main simulation
// ============================================================
int main(int argc, char* argv[])
{
  olb::initialize(&argc, &argv);
  OstreamManager clout(std::cout, "eel3dof");
  singleton::directories().setOutputDir("./tmp/");
  const std::string baseLogDir = singleton::directories().getLogOutDir();
  if (!eelEnsureDirectoryTree(baseLogDir)) {
    clout << "ERROR: could not create output directory: " << baseLogDir << std::endl;
    return 1;
  }

  CliParseResult cli = parseCommandLineDetailed(argc, argv, baseLogDir);
  if (cli.helpRequested && cli.diagnostics.errors.empty()) {
    clout << commandLineUsage(argc > 0 ? argv[0] : "11_lbm_eel_3dof");
    return 0;
  }

  CliDiagnostics configDiagnostics = validateRunConfig(cli.config);
  cli.diagnostics.errors.insert(cli.diagnostics.errors.end(),
                                configDiagnostics.errors.begin(),
                                configDiagnostics.errors.end());
  cli.diagnostics.warnings.insert(cli.diagnostics.warnings.end(),
                                  configDiagnostics.warnings.begin(),
                                  configDiagnostics.warnings.end());
  for (const std::string& warning : cli.diagnostics.warnings) {
    clout << "WARNING: " << warning << std::endl;
  }
  if (!cli.diagnostics.errors.empty()) {
    for (const std::string& error : cli.diagnostics.errors) {
      clout << "ERROR: " << error << std::endl;
    }
    clout << commandLineUsage(argc > 0 ? argv[0] : "11_lbm_eel_3dof");
    return 1;
  }

  RunConfig config = cli.config;
  EelParams& p = config.p;
  RunMode& runMode = config.runMode;
  StudyMode& studyMode = config.studyMode;
  SimulationCase& simCase = config.simCase;
  WarmupMode& warmupMode = config.warmupMode;
  BodyKinematics& bodyKinematics = config.bodyKinematics;
  bool& softBackboneDynamics = config.softBackboneDynamics;
  int& softBackboneCouplingIterations =
    config.softBackboneCouplingIterations;
  T& alphaIBM = config.alphaIBM;
  int& ibmIterations = config.ibmIterations;
  bool& rawGeometryOverride = config.rawGeometryOverride;
  bool& aspectGeometryOverride = config.aspectGeometryOverride;
  T& tCut = config.tCut;
  WallBoundary& wallBoundary = config.wallBoundary;
  std::string& summaryCsv = config.summaryCsv;
  std::string& sensitivityCsv = config.sensitivityCsv;
  std::string& runTag = config.runTag;
  bool& exportVelocity = config.exportVelocity;
  bool& exportVorticity = config.exportVorticity;
  bool& exportDiagnostics = config.exportDiagnostics;
  bool& exportBody = config.exportBody;

  if (rawGeometryOverride && !aspectGeometryOverride) {
    p.useAspectRatioGeometry = false;
  }
  if (!updateGeometryFromAspectRatio(p)) {
    clout << "ERROR: invalid aspect-ratio geometry. Require aspectRatio > 1 "
          << "and bodyAreaTarget > 0." << std::endl;
    return 1;
  }
  if (!p.useAspectRatioGeometry) {
    p.aspectRatio = p.totalGeometricLengthLU() / p.bodyWidthLU();
    p.bodyAreaTarget = p.capsuleAreaLU();
  }

  if (tCut < T(0)) {
    const T periodsForAverage = (p.eelFreq > T(1e-12)) ? T(5) / p.eelFreq : T(3);
    tCut = p.restTime + p.rampTime + std::max(T(3), periodsForAverage);
  }

  const std::string runsRootDir = baseLogDir + "runs/";
  if (!eelEnsureDirectoryTree(runsRootDir)) {
    clout << "ERROR: could not create run output root: " << runsRootDir << std::endl;
    return 1;
  }
  const std::string baseRunId =
    makeBaseRunId(p, simCase, config.bodyKinematics, runTag);
  const std::string runId = makeUniqueRunId(baseRunId, runsRootDir);
  const std::string runOutDir = runsRootDir + runId + "/";
  const std::string runVtkDir = runOutDir + "vtkData/";
  const std::string bodyOutDir = runOutDir + "bodyData/";
  const std::string velDir  = runVtkDir + "velocity/";
  const std::string vortDir = runVtkDir + "vorticity/";
  const std::string diagDir = runVtkDir + "diagnostics/";
  if (!eelEnsureDirectoryTree(runOutDir) ||
      (exportBody && !eelEnsureDirectoryTree(bodyOutDir)) ||
      (exportVelocity && !eelEnsureDirectoryTree(velDir)) ||
      (exportVorticity && !eelEnsureDirectoryTree(vortDir)) ||
      (exportDiagnostics && runMode == RunMode::Full &&
       !eelEnsureDirectoryTree(diagDir))) {
    clout << "ERROR: could not create per-run output directories under "
          << runOutDir << std::endl;
    return 1;
  }

  // Export cadence: preview writes every 5th frame, others every frame
  const int exportInterval = (runMode == RunMode::Preview) ? 5 : 1;
  const bool verificationMode = (studyMode == StudyMode::Verification);
  const bool writeSpatialExports = !verificationMode;

  const int nx = p.nx;
  const int ny = p.ny;
  const T dtLbm = p.dtLbm();
  const int nFrames = p.nFrames();

  clout << "==== Self-propelled eel 3-DOF (OpenLB-native fluid + custom IBM) ====" << std::endl;
  clout << "Run ID: " << runId << std::endl;
  clout << "Run output: " << runOutDir << std::endl;
  clout << "Study mode: " << studyModeName(studyMode) << std::endl;
  clout << "Sensitivity summary: " << sensitivityCsv << std::endl;
  clout << "Grid: " << nx << " x " << ny << "  tau=" << p.tau
        << "  nu=" << p.nu() << std::endl;
  clout << "Case: " << simulationCaseName(simCase)
        << "  wallBoundary=" << wallBoundaryName(wallBoundary) << std::endl;
  clout << "AR study: useAspectRatioGeometry=" << (p.useAspectRatioGeometry ? "true" : "false")
        << "  AR=" << p.aspectRatio
        << "  bodyAreaTarget=" << p.bodyAreaTarget
        << "  tCut=" << tCut << std::endl;
  clout << "Mode note: surge_only advances only the streamwise rigid-body DOF; "
        << "heave/yaw are excluded by model definition." << std::endl;
  clout << "Eel: centerlineL=" << p.centerlineLengthLU()
        << " lu  totalL=" << p.totalGeometricLengthLU()
        << " lu  W=" << p.bodyWidthLU()
        << " lu  freq=" << p.eelFreq
        << "  lambda=" << p.eelLambda
        << "  waveDirection=" << waveDirectionName(p.waveDirection)
        << "  bodyKinematics="
        << bodyKinematicsName(bodyKinematics) << std::endl;
  clout << "Warmup: mode=" << warmupModeName(warmupMode)
        << "  nWarmup=" << p.nWarmup << std::endl;
  clout << "Time: Ttotal=" << p.Ttotal << "  dt_anim=" << p.dtAnim
        << "  substeps=" << p.substeps << "  dt_lbm=" << dtLbm << std::endl;
  clout << "OpenMP threads: " << omp_get_max_threads()
        << "  mode=" << (runMode == RunMode::Preview ? "preview" :
                         runMode == RunMode::Standard ? "standard" : "full")
        << "  exportInterval=" << exportInterval << std::endl;
  clout << "Exports: velocity=" << (exportVelocity ? 1 : 0)
        << "  vorticity=" << (exportVorticity ? 1 : 0)
        << "  diagnostics=" << (exportDiagnostics ? 1 : 0)
        << " (full mode only)"
        << "  body=" << (exportBody ? 1 : 0) << std::endl;
  if (verificationMode) {
    clout << "Verification mode keeps the same solver physics/update order and "
          << "suppresses VTK/body exports to make parameter sweeps lighter."
          << std::endl;
  }
  clout << "IBM: alphaIBM=" << alphaIBM
        << "  ibmIterations=" << ibmIterations
        << " (forcing law: F = alphaIBM * rho_local * (Ud - U_interp) / dt;"
        << " rho_local=1.0, dt=1.0 in lattice units;"
        << " iterations use sparse Eulerian correction re-interpolation)"
        << std::endl;
  clout << "IBM legacy params: nIbmIters(requested)=" << config.legacyNIbmIters
        << "  spongeWidth=" << p.spongeWidth
        << "  spongeStrength=" << p.spongeStrength << std::endl;
  clout << "IBM warnings: meanSlip>" << p.warnMeanSlip
        << "  maxSlip>" << p.warnMaxSlip
        << "  maxMarkerForce>" << p.warnMarkerForce << std::endl;
  if (config.legacyNIbmIters != 1) {
    clout << "NOTE: --nIbmIters is retained for input compatibility, but "
          << "the corrected coupling uses one IBM force evaluation per fluid step."
          << std::endl;
  }

  // ---- OpenLB-native fluid lattice ----
  UnitConverter<T, DESCRIPTOR> converter(
    1.0,               // lattice units: dx = 1
    1.0,               // lattice units: dt = 1
    p.totalGeometricLengthLU(),  // characteristic length: capsule tip-to-tip length
    0.05,              // characteristic lattice velocity for reporting
    p.nu(),
    1.0
  );
  converter.print();

  IndicatorCuboid2D<T> domain({T(nx - 1), T(ny - 1)}, {T(0), T(0)});
  CuboidDecomposition2D<T> cuboidDecomposition(domain, 1.0, singleton::mpi().getSize());
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer, 3);

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1, 1});

  // Channel layout: top/bottom are no-slip walls; both ends are open
  // pressure boundaries for the self-propelled swimming case.
  IndicatorCuboid2D<T> leftOpen({0.25, T(ny - 1)}, {T(0), T(0)});
  IndicatorCuboid2D<T> rightOpen({0.25, T(ny - 1)}, {T(nx - 1) - 0.25, T(0)});
  superGeometry.rename(2, 3, 1, leftOpen);
  superGeometry.rename(2, 4, 1, rightOpen);
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();

  SuperLattice<T, DESCRIPTOR> sLattice(converter, superGeometry);
  dynamics::set<ForcedBGKdynamics<T, DESCRIPTOR>>(sLattice, superGeometry, 1);
  // Top/bottom walls (material 2):
  //   noslip   -> BounceBack (legacy default; channel walls)
  //   freeslip -> FullSlip   (specular reflection; removes channel-blockage
  //                            drag, important for AR shape comparisons)
  if (wallBoundary == WallBoundary::FreeSlip) {
    boundary::set<boundary::FullSlip>(sLattice, superGeometry, 2);
  } else {
    boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  }
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
  sLattice.template setParameter<descriptors::OMEGA>(p.omega_lbm());

  AnalyticalConst2D<T,T> rho0(1.0);
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> f0(0.0, 0.0);
  auto fluidAndBoundary = superGeometry.getMaterialIndicator({1, 2, 3, 4});
  sLattice.iniEquilibrium(fluidAndBoundary, rho0, u0);
  sLattice.defineRhoU(fluidAndBoundary, rho0, u0);
  fields::set<descriptors::FORCE>(sLattice, fluidAndBoundary, f0);
  sLattice.initialize();

  // ---- Sponge zone near open (left/right) boundaries ----
  //
  // The previous solver blended raw populations toward f_eq.  In the OpenLB
  // workflow we keep this as a lattice FORCE-field damping layer, so the
  // damping is applied by ForcedBGKdynamics instead of direct population edits.
  auto resetEulerForce = [&]() {
    resetEulerForceField(p, sLattice, /*fixedInflow=*/false);
  };
  resetEulerForce();

  // ============================================================
  //  3-DOF rigid-body state  (all in LBM natural units, δt = 1)
  //
  //  Velocities : Vx, Vy   [lu / step]
  //  Angular vel: omegaZ    [rad / step]
  //  Positions  : xCm, yCm  [lu]
  //  Orientation: theta      [rad]
  //
  //  Newton-Euler per step:
  //    ΔV = F_body / M       — no extra dt factor (δt = 1)
  //    Δx = V                — no extra dt factor
  //    Δθ = ω                — no extra dt factor
  // ============================================================
  // Initial placement is a setup parameter only; it does not alter the
  // swimmer kinematics or coupling.  The right-margin clamp below keeps
  // the body far enough from the outlet for wake visualization.
  BodyState bodyState;
  bodyState.xCm = T(nx) * p.initialPositionFactor;
  bodyState.yCm = T(ny) * 0.5;
  bodyState.theta = 0.0;
  bodyState.Vx = 0.0;
  bodyState.Vy = 0.0;
  bodyState.omegaZ = 0.0;
  T& xCm = bodyState.xCm;
  T& yCm = bodyState.yCm;
  T& theta = bodyState.theta;
  T& Vx = bodyState.Vx;
  T& Vy = bodyState.Vy;
  T& omegaZ = bodyState.omegaZ;

  // Body inertial properties.  The centerline length is used by the gait
  // construction; the capsule tip-to-tip length is the characteristic body
  // length for mass/inertia/Re diagnostics.  Density now comes from the
  // selected body material (Dragon Skin 20 by default).
  const MaterialProperties material = resolveMaterialProperties(p);
  const PlanarRodSectionEstimate rodSection =
    estimatePlanarRodSection(p, material);
  const SoftBackboneConfig softBackbone =
    makeSoftBackboneConfig(p, rodSection, std::max(2, p.nSpine - 1),
                           config.softBackboneAddedMassFrac);
  BodyInertia bodyInertia = computeBodyInertia(p, material.densityRatio);
  T bodyCenterlineLength = bodyInertia.bodyCenterlineLength;
  T bodyLength = bodyInertia.bodyLength;
  T bodyWidth  = bodyInertia.bodyWidth;
  T bodyArea   = bodyInertia.bodyArea;
  T mass       = bodyInertia.mass;
  T Ibody      = bodyInertia.Ibody;
  T mAddedSurge_theory = bodyInertia.mAddedSurgeTheory;
  T mAddedHeave_theory = bodyInertia.mAddedHeaveTheory;
  T Iadded_theory = bodyInertia.IaddedTheory;
  T Msurge      = bodyInertia.Msurge;
  T Mheave      = bodyInertia.Mheave;
  T Itotal      = bodyInertia.Itotal;

  // ClampAudit separates a one-shot setup adjustment from clamps that fire
  // during the warmup or main loop.  initialPlacement* tracks pre-warmup
  // adjustments to xCm caused by initialPositionFactor / rightMargin.
  // runtimeDomainClamp* tracks clampBodyBeforeMarkers() corrections that
  // happen once the simulation is actually advancing.  This avoids the
  // misleading "domainClampHit=1, domainClampCount=1" reading produced by a
  // benign one-time setup nudge.
  ClampAudit clampAudit;

  const T rightMargin = T(5) * bodyLength;
  const T unclampedXCm = xCm;
  const T maxInitialXCm = T(nx) - rightMargin;
  if (maxInitialXCm <= T(0)) {
    clout << "ERROR: domain is too short for rightMargin=5*bodyLength. "
          << "nx=" << nx << "  bodyLength=" << bodyLength
          << "  rightMargin=" << rightMargin << std::endl;
    return 1;
  }
  xCm = std::min(xCm, maxInitialXCm);
  if (xCm != unclampedXCm) {
    clampAudit.initialPlacementClamped = true;
    ++clampAudit.initialPlacementClampCount;
  }
  T xCmRef    = xCm;
  T yCmRef    = yCm;
  T thetaRef  = theta;

  clout << "Initial placement: factor=" << p.initialPositionFactor
        << "  xCm=" << xCm
        << "  nx=" << nx
        << "  distanceToRight=" << (T(nx) - xCm)
        << "  rightMargin=" << rightMargin;
  if (xCm != unclampedXCm) {
    clout << "  (clamped from " << unclampedXCm << ")";
  }
  clout << std::endl;

  // ----------------------------------------------------------------
  //  Added mass — 2D potential-flow theory for an ellipse, scaled by
  //  addedMassFrac.  With addedMassFrac = 0 (default) the IBM fluid
  //  forces handle all inertial coupling; no explicit added mass.
  // ----------------------------------------------------------------
  clout << "Body: centerlineL=" << bodyCenterlineLength
        << "  totalL=" << bodyLength << "  W=" << bodyWidth
        << "  capsuleArea=" << bodyArea << std::endl;
  clout << "Material: " << material.name
        << "  rhoBodyRatio=" << material.densityRatio
        << "  rhoBodySI=" << material.densityKgM3 << " kg/m^3"
        << "  E=" << material.youngModulusPa << " Pa"
        << "  nu=" << material.poissonRatio
        << "  G=" << material.shearModulusPa << " Pa"
        << "  K=" << material.bulkModulusPa << " Pa"
        << "  dampingRatio=" << material.dampingRatio << std::endl;
  if (rodSection.valid) {
    clout << "Soft rod estimate: physicalL=" << rodSection.physicalLengthM
          << " m  width=" << rodSection.widthM
          << " m  thickness=" << rodSection.thicknessM
          << " m  EI=" << rodSection.bendingStiffnessNm2
          << " N*m^2  EA=" << rodSection.axialStiffnessN
          << " N  m/L=" << rodSection.massPerLengthKgM
          << " kg/m" << std::endl;
  } else {
    clout << "Soft rod estimate: inactive. Set --physicalBodyLengthM=<m> "
          << "to map Dragon Skin 20 stiffness to the lattice geometry."
          << std::endl;
  }
  if (softBackbone.valid) {
    const T probePeriod =
      (p.eelFreq > T(1e-12)) ? T(0.25) / p.eelFreq : T(0);
    const T probeTime = p.restTime + p.rampTime + probePeriod;
    const SoftBackboneCurvatureProfile preferredBackbone =
      preferredBackboneCurvatureWave(p, softBackbone, probeTime, p.dtAnim);
    const SoftBackboneState straightBackbone =
      makeStraightBackboneState(softBackbone.nSegments);
    const SoftBackboneTorqueResult straightBackboneTorque =
      computeSoftBackboneTorques(softBackbone, straightBackbone,
                                 preferredBackbone);
    const T addedRatio =
      (softBackbone.segmentRotationalInertiaKgM2 > T(0))
        ? softBackbone.addedSegmentRotationalInertiaKgM2 /
          softBackbone.segmentRotationalInertiaKgM2
        : T(0);
    clout << "Soft backbone foundation: segments=" << softBackbone.nSegments
          << "  centerlineL=" << softBackbone.centerlineLengthM
          << " m  ds=" << softBackbone.dsM
          << " m  Ktheta=" << softBackbone.jointAngleStiffnessNm
          << " N*m/rad  Ctheta=" << softBackbone.jointAngleDampingNms
          << " N*m*s/rad  Iseg=" << softBackbone.segmentRotationalInertiaKgM2
          << " kg*m^2  IsegAdded=" << softBackbone.addedSegmentRotationalInertiaKgM2
          << " kg*m^2  (frac="
          << config.softBackboneAddedMassFrac
          << ", I_added/I_struct=" << addedRatio
          << ")  maxPreferredCurvature="
          << preferredBackbone.maxAbsCurvature
          << " 1/m  maxStraightMoment="
          << straightBackboneTorque.maxAbsJointMomentNm
          << " N*m" << std::endl;
  } else {
    clout << "Soft backbone foundation: inactive; rod estimate is invalid."
          << std::endl;
  }
  if (bodyKinematics == BodyKinematics::SoftBackbone &&
      !softBackbone.valid) {
    clout << "ERROR: --bodyKinematics=soft_backbone needs a valid "
          << "Dragon Skin soft-rod estimate. Check --physicalBodyLengthM, "
          << "--bodyThicknessM, and material properties." << std::endl;
    return 1;
  }
  if (bodyKinematics == BodyKinematics::SoftBackbone) {
    clout << "Soft backbone coupling: dynamics="
          << (softBackboneDynamics ? "enabled" : "disabled")
          << "  integrator=implicit_mean_orientation_affine_load_filter"
          << "  fluidTorqueScale="
          << config.softBackboneFluidTorqueScale
          << "  fluidTorqueFilterTime="
          << config.softBackboneFluidTorqueFilterTime
          << "  maxAngleStep="
          << config.softBackboneMaxAngleStep
          << " rad/step"
          << "  couplingIterations="
          << softBackboneCouplingIterations
          << "  couplingRelaxation="
          << config.softBackboneCouplingRelaxation
          << "  couplingTolerance="
          << config.softBackboneCouplingTolerance
          << "  loadProjection="
          << softBackboneLoadProjectionName(
               config.softBackboneLoadProjection);
    if (softBackboneDynamics) {
      clout << "  abortOnInstability="
            << (config.softBackboneAbortOnInstability ? 1 : 0)
            << "  abortMeanSlip="
            << config.softBackboneAbortMeanSlip
            << "  abortMaxSlip="
            << config.softBackboneAbortMaxSlip
            << "  abortSaturatedFrames="
            << config.softBackboneAbortSaturatedFrames;
    }
    clout << std::endl;
  }
  clout << "mass=" << mass << "  Ibody=" << Ibody << std::endl;
  clout << "addedMassFrac=" << p.addedMassFrac
        << "  mAdded_surge_theory=" << mAddedSurge_theory
        << "  mAdded_heave_theory=" << mAddedHeave_theory
        << "  Iadded_theory=" << Iadded_theory << std::endl;
  clout << "Msurge=" << Msurge << "  Mheave=" << Mheave
        << "  Itotal=" << Itotal << std::endl;

  // ----------------------------------------------------------------
  //  Non-dimensional / resolution diagnostics
  // ----------------------------------------------------------------
  T tailAmpLU = p.tailAmplitudeLU();
  T tailPP    = p.tailPeakToPeakLU();
  T freqLat   = p.eelFreq * dtLbm;  // frequency in [1/step]
  T blockageRatio = bodyWidth / T(ny);
  T clearanceLU   = 0.5 * (T(ny) - bodyWidth);

  clout << "---- resolution diagnostics ----" << std::endl;
  clout << "  body length / dx     = " << bodyLength << "  (aim >= 100)" << std::endl;
  clout << "  body width  / dx     = " << bodyWidth  << "  (aim >= 10)" << std::endl;
  clout << "  tail amplitude       = " << tailAmpLU  << " lu  ("
        << (tailAmpLU / bodyCenterlineLength) << " centerline L)" << std::endl;
  clout << "  tail peak-to-peak    = " << tailPP << " lu" << std::endl;
  clout << "  blockage ratio W/ny  = " << blockageRatio
        << "  (aim < 0.05)" << std::endl;
  clout << "  lateral clearance    = " << clearanceLU << " lu  ("
        << (clearanceLU / bodyLength) << " L)" << std::endl;
  clout << "  freq_lattice         = " << freqLat << " /step" << std::endl;
  clout << "  steps/undulation     = " << (1.0 / freqLat) << std::endl;

  // surge_only: reduced-order self-propelled model — only streamwise
  // translation is a dynamic state.  Heave/yaw are pinned by definition.
  auto deformationAmpScale = []() -> T { return T(1); };

  SoftBackboneDynamicsParams softDynamicsParams;
  softDynamicsParams.fluidTorqueScale = config.softBackboneFluidTorqueScale;
  softDynamicsParams.maxAngleStep = config.softBackboneMaxAngleStep;
  SoftBackboneState softDynamicState;
  std::vector<T> filteredSoftSegmentTorqueNm;
  bool softTorqueFilterInitialized = false;
  if (bodyKinematics == BodyKinematics::SoftBackbone) {
    softDynamicState =
      preferredBackboneStateWave(p, softBackbone, T(0), dtLbm,
                                 deformationAmpScale());
  }

  auto applyModeDefinitionState = [&]() {
    enforceModeDefinitionState({xCmRef, yCmRef, thetaRef}, bodyState);
  };

  bool runtimeDomainClampWarned = false;
  const T markerDomainGuard = T(2);
  const T bodyClampSafety = T(2);
  auto clampBodyBeforeMarkers = [&]() -> bool {
    const T halfLength = 0.5 * bodyLength;
    const T halfWidthWithWave = p.bodyRadius + tailAmpLU
                              + markerDomainGuard + bodyClampSafety;
    const T c = std::abs(std::cos(theta));
    const T s = std::abs(std::sin(theta));
    const T xMargin = c * halfLength + s * halfWidthWithWave;
    const T yMargin = s * halfLength + c * halfWidthWithWave;
    if (T(nx - 1) <= 2.0 * xMargin ||
        T(ny - 1) <= 2.0 * yMargin) {
      clout << "ERROR: domain is too small for the capsule plus IBM support. "
            << "xMargin=" << xMargin << " yMargin=" << yMargin
            << " nx=" << nx << " ny=" << ny << std::endl;
      return false;
    }

    const T oldX = xCm;
    const T oldY = yCm;
    xCm = std::clamp(xCm, xMargin, T(nx - 1) - xMargin);
    yCm = std::clamp(yCm, yMargin, T(ny - 1) - yMargin);

    if (xCm != oldX) Vx = T(0);
    if (xCm != oldX || yCm != oldY) {
      clampAudit.runtimeDomainClampHit = true;
      ++clampAudit.runtimeDomainClampCount;
      if (!runtimeDomainClampWarned) {
        clout << "WARNING: clamped body CM before marker generation from ("
              << oldX << ", " << oldY << ") to ("
              << xCm << ", " << yCm << ")." << std::endl;
        runtimeDomainClampWarned = true;
      }
      if (yCm != oldY) {
        yCmRef = yCm;
      }
    }
    return true;
  };

  auto markersInDomain = [&](const LagrangianMarkers& mk, const char* context) -> bool {
    if (mk.size() == 0) return false;
    T minX = mk.x[0], maxX = mk.x[0], minY = mk.y[0], maxY = mk.y[0];
    for (int m = 1; m < mk.size(); ++m) {
      minX = std::min(minX, mk.x[m]);
      maxX = std::max(maxX, mk.x[m]);
      minY = std::min(minY, mk.y[m]);
      maxY = std::max(maxY, mk.y[m]);
    }
    const bool ok = (minX >= markerDomainGuard &&
                     maxX <= T(nx - 1) - markerDomainGuard &&
                     minY >= markerDomainGuard &&
                     maxY <= T(ny - 1) - markerDomainGuard);
    if (!ok) {
      clout << "ERROR: IBM markers left the lattice during " << context
            << ". bounds x=[" << minX << ", " << maxX << "] y=["
            << minY << ", " << maxY << "] domain=[0," << (nx - 1)
            << "]x[0," << (ny - 1) << "]. Terminating cleanly."
            << std::endl;
    }
    return ok;
  };

  auto buildMarkersFromSoftStateChecked =
      [&](const SoftBackboneState& state, LagrangianMarkers& target,
          const char* context) -> bool {
    applyModeDefinitionState();
    if (!clampBodyBeforeMarkers()) return false;
    applyModeDefinitionState();
    buildLagrangianMarkersFromSoftBackboneState(
      p, softBackbone, state, Vx, Vy, omegaZ,
      xCm, yCm, theta, dtLbm, target);
    return markersInDomain(target, context);
  };

  auto buildMarkersChecked = [&](T t, LagrangianMarkers& target,
                                 const char* context) -> bool {
    if (bodyKinematics == BodyKinematics::SoftBackbone &&
        softBackboneDynamics) {
      return buildMarkersFromSoftStateChecked(
        softDynamicState, target, context);
    }
    applyModeDefinitionState();
    if (!clampBodyBeforeMarkers()) return false;
    applyModeDefinitionState();
    if (bodyKinematics == BodyKinematics::SoftBackbone) {
      buildLagrangianMarkersFromSoftBackbone(
        p, softBackbone, t, Vx, Vy, omegaZ, xCm, yCm, theta, dtLbm,
        deformationAmpScale(), target);
    } else {
      buildLagrangianMarkers(
        p, t, Vx, Vy, omegaZ, xCm, yCm, theta, dtLbm,
        deformationAmpScale(), target);
    }
    return markersInDomain(target, context);
  };

  // ---- Build initial markers ----
  LagrangianMarkers markers;
  if (!buildMarkersChecked(0.0, markers, "initialization")) {
    return 1;
  }

  // Kinematic time offset: gait clock starts at 0 at history t=0 so per-cycle
  // averages from cycle 1 are phase-aligned with the gait kinematics.
  auto kinematicTimeAt = [&](T reportedTime) -> T {
    return reportedTime;
  };

  // ---- Warmup ----
  IbmResult ibmRes;
  if (warmupMode == WarmupMode::None) {
    clout << "Warmup mode: none — skipping warmup loop." << std::endl;
  } else {
    clout << "Starting warmup (" << p.nWarmup << " steps, mode=rest)..."
          << std::endl;
    for (int step = 0; step < p.nWarmup; ++step) {
      // Rest warmup: actuationProfile at tWarmup=0 returns ampScale=0, so the
      // body stays in its still pose and the fluid relaxes to a quiescent /
      // boundary-consistent state.
      if (!buildMarkersChecked(T(0), markers, "warmup")) {
        return 1;
      }
      ibmStep(p, markers, xCm, yCm, Vx, Vy, omegaZ, sLattice,
              alphaIBM, ibmIterations, resetEulerForce, ibmRes);
      sLattice.collideAndStream();
      resetEulerForce();

      if (step % 100 == 0) {
        clout << "  warmup " << step << "/" << p.nWarmup << std::endl;
      }
    }
    clout << "Warmup complete." << std::endl;
  }

  // ============================================================
  //  History arrays
  // ============================================================
  SimulationHistory history;

  // ============================================================
  //  Output helpers
  // ============================================================

  ImageDataVtiWriter vtiWriter(nx, ny);
  ExportSnapshot exportFields(nx, ny);

  // ============================================================
  //  Main simulation loop
  // ============================================================
  clout << "Starting main simulation: " << nFrames << " frames, "
        << p.substeps << " substeps each." << std::endl;

  int softBackboneSaturatedFrameCount = 0;

  for (int frame = 0; frame < nFrames; ++frame) {
    T tStart = frame * p.dtAnim;

    T FxAccum = 0.0, FyAccum = 0.0, TzAccum = 0.0, PAccum = 0.0;
    // v6 per-frame accumulators for the power decomposition and the
    // marker-local residual slip from MDF iterations.
    T PRigidAccum = 0.0, PDefAccum = 0.0;
    T meanSlipAccum = 0.0;
    T meanResidualSlipAccum = 0.0;
    T meanMarkerForceAccum = 0.0;
    T meanDesiredSpeedAccum = 0.0;
    T maxSlipFrame = 0.0;
    T maxResidualSlipFrame = 0.0;
    T maxMarkerForceFrame = 0.0;
    T softFluidTorqueAccum = 0.0;
    T softAngleStepAccum = 0.0;
    T maxSoftFluidTorqueFrame = 0.0;
    T maxSoftAngleStepFrame = 0.0;
    T softCouplingResidualAccum = 0.0;
    T maxSoftCouplingResidualFrame = 0.0;
    T softCouplingItersAccum = 0.0;
    T maxSoftCouplingItersFrame = 0.0;
    T softElasticEnergyAccum = 0.0;
    T maxSoftElasticEnergyFrame = 0.0;
    T softDampingPowerAccum = 0.0;
    T softActuatorPowerAccum = 0.0;
    T softAbsActuatorPowerAccum = 0.0;
    T maxAbsSoftActuatorPowerFrame = 0.0;
    T softAppliedFluidPowerAccum = 0.0;
    T softAbsAppliedFluidPowerAccum = 0.0;
    T maxAbsSoftAppliedFluidPowerFrame = 0.0;
    int softDiagSamples = 0;
    int softCouplingSamples = 0;
    int ibmDiagSamples = 0;

    auto accumulateIbmDiagnostics = [&]() {
      if (std::isfinite(ibmRes.meanSlipMag) &&
          std::isfinite(ibmRes.maxSlipMag) &&
          std::isfinite(ibmRes.meanMarkerForceMag) &&
          std::isfinite(ibmRes.maxMarkerForceMag) &&
          std::isfinite(ibmRes.meanDesiredMarkerSpeedMag)) {
        meanSlipAccum += ibmRes.meanSlipMag;
        meanMarkerForceAccum += ibmRes.meanMarkerForceMag;
        meanDesiredSpeedAccum += ibmRes.meanDesiredMarkerSpeedMag;
        maxSlipFrame = std::max(maxSlipFrame, ibmRes.maxSlipMag);
        maxMarkerForceFrame =
          std::max(maxMarkerForceFrame, ibmRes.maxMarkerForceMag);
        if (std::isfinite(ibmRes.meanResidualSlipMag) &&
            std::isfinite(ibmRes.maxResidualSlipMag)) {
          meanResidualSlipAccum += ibmRes.meanResidualSlipMag;
          maxResidualSlipFrame =
            std::max(maxResidualSlipFrame, ibmRes.maxResidualSlipMag);
        }
        ++ibmDiagSamples;
      }
    };

    auto accumulateSoftDiagnostics =
        [&](const SoftBackboneDynamicsDiagnostics& softDynDiag) {
      if (std::isfinite(softDynDiag.maxAbsFluidSegmentTorqueNm) &&
          std::isfinite(softDynDiag.maxAbsAngleStep) &&
          std::isfinite(softDynDiag.elasticEnergyJ) &&
          std::isfinite(softDynDiag.dampingPowerW) &&
          std::isfinite(softDynDiag.actuatorPowerProxyW) &&
          std::isfinite(softDynDiag.absActuatorPowerProxyW) &&
          std::isfinite(softDynDiag.appliedFluidPowerW) &&
          std::isfinite(softDynDiag.absAppliedFluidPowerW)) {
        softFluidTorqueAccum += softDynDiag.maxAbsFluidSegmentTorqueNm;
        softAngleStepAccum += softDynDiag.maxAbsAngleStep;
        maxSoftFluidTorqueFrame =
          std::max(maxSoftFluidTorqueFrame,
                   softDynDiag.maxAbsFluidSegmentTorqueNm);
        maxSoftAngleStepFrame =
          std::max(maxSoftAngleStepFrame,
                   softDynDiag.maxAbsAngleStep);
        softElasticEnergyAccum += softDynDiag.elasticEnergyJ;
        maxSoftElasticEnergyFrame =
          std::max(maxSoftElasticEnergyFrame,
                   softDynDiag.elasticEnergyJ);
        softDampingPowerAccum += softDynDiag.dampingPowerW;
        softActuatorPowerAccum += softDynDiag.actuatorPowerProxyW;
        softAbsActuatorPowerAccum += softDynDiag.absActuatorPowerProxyW;
        maxAbsSoftActuatorPowerFrame =
          std::max(maxAbsSoftActuatorPowerFrame,
                   softDynDiag.absActuatorPowerProxyW);
        softAppliedFluidPowerAccum += softDynDiag.appliedFluidPowerW;
        softAbsAppliedFluidPowerAccum += softDynDiag.absAppliedFluidPowerW;
        maxAbsSoftAppliedFluidPowerFrame =
          std::max(maxAbsSoftAppliedFluidPowerFrame,
                   softDynDiag.absAppliedFluidPowerW);
        ++softDiagSamples;
      }
    };

    auto wrapAngleLocal = [](T angle) -> T {
      angle = std::fmod(angle + M_PI, T(2) * M_PI);
      if (angle < T(0)) angle += T(2) * M_PI;
      return angle - M_PI;
    };

    auto backboneStateResidual =
        [&](const SoftBackboneState& previous,
            const SoftBackboneState& next) -> T {
      if (previous.theta.size() != next.theta.size() ||
          previous.omega.size() != next.omega.size()) {
        return std::numeric_limits<T>::infinity();
      }
      T residual = T(0);
      for (size_t i = 0; i < previous.theta.size(); ++i) {
        const T dTheta = std::abs(
          wrapAngleLocal(next.theta[i] - previous.theta[i]));
        const T dOmegaAsAngle =
          std::abs(next.omega[i] - previous.omega[i]) * dtLbm;
        residual = std::max(residual, std::max(dTheta, dOmegaAsAngle));
      }
      return residual;
    };

    auto relaxedBackboneState =
        [&](const SoftBackboneState& previous,
            const SoftBackboneState& next) -> SoftBackboneState {
      const T relax = config.softBackboneCouplingRelaxation;
      if (relax >= T(1)) {
        return next;
      }
      SoftBackboneState out = previous;
      if (out.theta.size() != next.theta.size() ||
          out.omega.size() != next.omega.size()) {
        return next;
      }
      for (size_t i = 0; i < out.theta.size(); ++i) {
        out.theta[i] =
          wrapAngleLocal(out.theta[i] +
                         relax * wrapAngleLocal(next.theta[i] - out.theta[i]));
        out.omega[i] += relax * (next.omega[i] - out.omega[i]);
      }
      return out;
    };

    auto accumulateSoftCouplingDiagnostics =
        [&](T residual, int itersUsed) {
      if (itersUsed > 0 && std::isfinite(residual)) {
        softCouplingResidualAccum += residual;
        maxSoftCouplingResidualFrame =
          std::max(maxSoftCouplingResidualFrame, residual);
        softCouplingItersAccum += T(itersUsed);
        maxSoftCouplingItersFrame =
          std::max(maxSoftCouplingItersFrame, T(itersUsed));
        ++softCouplingSamples;
      }
    };

    auto filteredSoftTorqueCandidate =
        [&](const std::vector<T>& rawTorque,
            const std::vector<T>& filterStart,
            bool filterInitializedAtStepStart) -> std::vector<T> {
      const T filterTime = config.softBackboneFluidTorqueFilterTime;
      if (!(filterTime > T(0)) || rawTorque.empty() ||
          !filterInitializedAtStepStart ||
          filterStart.size() != rawTorque.size()) {
        return rawTorque;
      }
      const T alpha = T(1) - std::exp(-dtLbm / filterTime);
      std::vector<T> filtered = rawTorque;
      for (size_t i = 0; i < rawTorque.size(); ++i) {
        filtered[i] = filterStart[i] + alpha * (rawTorque[i] - filterStart[i]);
      }
      return filtered;
    };

    auto commitFilteredSoftTorque =
        [&](const std::vector<T>& torqueForDynamics) {
      if (config.softBackboneFluidTorqueFilterTime > T(0) &&
          !torqueForDynamics.empty()) {
        filteredSoftSegmentTorqueNm = torqueForDynamics;
        softTorqueFilterInitialized = true;
      }
    };

    for (int sub = 0; sub < p.substeps; ++sub) {
      T tSub = tStart + sub * dtLbm;
      T tKinematics = kinematicTimeAt(tSub);

      // Coupling order: build geometry, evaluate IBM, optionally solve a
      // fixed-point soft-backbone update from the same t^n state, then stream
      // the final IBM force field once.
      if (softBackboneDynamics &&
          bodyKinematics == BodyKinematics::SoftBackbone &&
          softBackboneCouplingIterations > 1) {
        const SoftBackboneState stateAtStepStart = softDynamicState;
        const SoftBackboneState preferredNext =
          preferredBackboneStateWave(
            p, softBackbone, tKinematics + dtLbm, dtLbm,
            deformationAmpScale());
        SoftBackboneState guessState = stateAtStepStart;
        const std::vector<T> filterStart = filteredSoftSegmentTorqueNm;
        const bool filterInitializedAtStepStart = softTorqueFilterInitialized;
        std::vector<T> torqueForDynamics;
        SoftBackboneDynamicsDiagnostics softDynDiag;
        T couplingResidual = std::numeric_limits<T>::infinity();
        int couplingItersUsed = 0;
        for (int couplingIter = 0;
             couplingIter < softBackboneCouplingIterations;
             ++couplingIter) {
          const SoftBackboneState previousGuess = guessState;
          if (!buildMarkersFromSoftStateChecked(
                guessState, markers, "soft-backbone coupling")) {
            return 1;
          }
          ibmStep(p, markers, xCm, yCm, Vx, Vy, omegaZ, sLattice,
                  alphaIBM, ibmIterations, resetEulerForce, ibmRes);
          const SoftBackboneForceProjection softProjection =
            projectMarkerForcesToSoftBackbone(
              p, softBackbone, guessState, xCm, yCm, theta, markers,
              config.softBackboneLoadProjection);
          torqueForDynamics =
            filteredSoftTorqueCandidate(softProjection.segmentTorqueNm,
                                        filterStart,
                                        filterInitializedAtStepStart);
          SoftBackboneState nextState = stateAtStepStart;
          softDynDiag = advanceSoftBackboneImplicit(
            softBackbone, preferredNext, torqueForDynamics,
            dtLbm, softDynamicsParams, nextState);
          guessState = relaxedBackboneState(previousGuess, nextState);
          couplingResidual =
            backboneStateResidual(previousGuess, guessState);
          couplingItersUsed = couplingIter + 1;
          if (config.softBackboneCouplingTolerance > T(0) &&
              couplingResidual <= config.softBackboneCouplingTolerance) {
            break;
          }
        }
        commitFilteredSoftTorque(torqueForDynamics);
        softDynamicState = guessState;
        if (!buildMarkersFromSoftStateChecked(
              softDynamicState, markers, "soft-backbone coupling final")) {
          return 1;
        }
        ibmStep(p, markers, xCm, yCm, Vx, Vy, omegaZ, sLattice,
                alphaIBM, ibmIterations, resetEulerForce, ibmRes);
        accumulateIbmDiagnostics();
        accumulateSoftDiagnostics(softDynDiag);
        accumulateSoftCouplingDiagnostics(couplingResidual,
                                          couplingItersUsed);
      } else {
        if (!buildMarkersChecked(tKinematics, markers, "main loop")) {
          return 1;
        }
        ibmStep(p, markers, xCm, yCm, Vx, Vy, omegaZ, sLattice,
                alphaIBM, ibmIterations, resetEulerForce, ibmRes);
        accumulateIbmDiagnostics();
        if (softBackboneDynamics &&
            bodyKinematics == BodyKinematics::SoftBackbone) {
          const SoftBackboneForceProjection softProjection =
            projectMarkerForcesToSoftBackbone(
              p, softBackbone, softDynamicState, xCm, yCm, theta, markers,
              config.softBackboneLoadProjection);
          const std::vector<T> filterStart = filteredSoftSegmentTorqueNm;
          const bool filterInitializedAtStepStart =
            softTorqueFilterInitialized;
          const std::vector<T> torqueForDynamics =
            filteredSoftTorqueCandidate(softProjection.segmentTorqueNm,
                                        filterStart,
                                        filterInitializedAtStepStart);
          commitFilteredSoftTorque(torqueForDynamics);
          const SoftBackboneState preferredNext =
            preferredBackboneStateWave(
              p, softBackbone, tKinematics + dtLbm, dtLbm,
              deformationAmpScale());
          const SoftBackboneDynamicsDiagnostics softDynDiag =
            advanceSoftBackboneImplicit(
              softBackbone, preferredNext, torqueForDynamics,
              dtLbm, softDynamicsParams, softDynamicState);
          accumulateSoftDiagnostics(softDynDiag);
          accumulateSoftCouplingDiagnostics(T(0), 1);
        }
      }
      sLattice.collideAndStream();
      // Clear the one-step IBM force after it has been consumed by OpenLB's
      // ForcedBGK collision.  The next IBM pass rebuilds the body force from
      // the updated OpenLB populations.
      resetEulerForce();

      T fxS = ibmRes.forceX;
      T fyS = ibmRes.forceY;
      T tzS = ibmRes.torque;
      T pS  = ibmRes.power;
      T pRigidS = ibmRes.powerRigid;
      T pDefS   = ibmRes.powerDef;

      bool valid = std::isfinite(fxS) && std::isfinite(fyS)
                && std::isfinite(tzS) && std::abs(fxS) < 500.0;
      if (valid) {
        FxAccum += fxS;
        FyAccum += fyS;
        TzAccum += tzS;
        PAccum  += pS;
        if (std::isfinite(pRigidS)) PRigidAccum += pRigidS;
        if (std::isfinite(pDefS))   PDefAccum   += pDefS;

        // Newton-Euler: ΔV = F/M  (δt = 1 in lattice units).
        // Forces are the body-reaction from this fluid step's IBM evaluation.
        // surge_only: only the streamwise force component is integrated.
        applyRigidBodyForceUpdate(tKinematics, p.restTime,
                                  fxS, bodyInertia, bodyState);
      }

      // Position/orientation: Δx = V  (δt = 1).  surge_only pins heave/yaw.
      advanceRigidBodyPose({xCmRef, yCmRef, thetaRef}, bodyState);

      if (!clampBodyBeforeMarkers()) {
        return 1;
      }
      applyModeDefinitionState();
    }

    // Frame averages
    T FxAvg = FxAccum / p.substeps;
    T FyAvg = FyAccum / p.substeps;
    T TzAvg = TzAccum / p.substeps;
    T PAvg  = PAccum  / p.substeps;
    T PRigidAvg = PRigidAccum / p.substeps;
    T PDefAvg   = PDefAccum   / p.substeps;
    T meanSlipFrame = 0.0;
    T meanResidualSlipFrame = 0.0;
    T meanMarkerForceFrame = 0.0;
    T meanDesiredSpeedFrame = 0.0;
    T normalizedSlipFrame = 0.0;
    T meanSoftFluidTorqueFrame = 0.0;
    T meanSoftAngleStepFrame = 0.0;
    T softCouplingResidualFrame = 0.0;
    T softCouplingItersUsedFrame = 0.0;
    T meanSoftElasticEnergyFrame = 0.0;
    T meanSoftDampingPowerFrame = 0.0;
    T meanSoftActuatorPowerProxyFrame = 0.0;
    T meanAbsSoftActuatorPowerProxyFrame = 0.0;
    T meanSoftAppliedFluidPowerFrame = 0.0;
    T meanAbsSoftAppliedFluidPowerFrame = 0.0;
    if (ibmDiagSamples > 0) {
      const T invDiag = T(1) / T(ibmDiagSamples);
      meanSlipFrame = meanSlipAccum * invDiag;
      meanResidualSlipFrame = meanResidualSlipAccum * invDiag;
      meanMarkerForceFrame = meanMarkerForceAccum * invDiag;
      meanDesiredSpeedFrame = meanDesiredSpeedAccum * invDiag;
      if (meanDesiredSpeedFrame > T(1e-12)) {
        normalizedSlipFrame = meanSlipFrame / meanDesiredSpeedFrame;
      }
    }
    if (softDiagSamples > 0) {
      const T invSoft = T(1) / T(softDiagSamples);
      meanSoftFluidTorqueFrame = softFluidTorqueAccum * invSoft;
      meanSoftAngleStepFrame = softAngleStepAccum * invSoft;
      meanSoftElasticEnergyFrame = softElasticEnergyAccum * invSoft;
      meanSoftDampingPowerFrame = softDampingPowerAccum * invSoft;
      meanSoftActuatorPowerProxyFrame =
        softActuatorPowerAccum * invSoft;
      meanAbsSoftActuatorPowerProxyFrame =
        softAbsActuatorPowerAccum * invSoft;
      meanSoftAppliedFluidPowerFrame =
        softAppliedFluidPowerAccum * invSoft;
      meanAbsSoftAppliedFluidPowerFrame =
        softAbsAppliedFluidPowerAccum * invSoft;
    }
    if (softCouplingSamples > 0) {
      const T invSoftCoupling = T(1) / T(softCouplingSamples);
      softCouplingResidualFrame =
        softCouplingResidualAccum * invSoftCoupling;
      softCouplingItersUsedFrame =
        softCouplingItersAccum * invSoftCoupling;
    }
    if (!std::isfinite(FxAvg) || !std::isfinite(FyAvg) ||
        !std::isfinite(TzAvg) || !std::isfinite(PAvg)) {
      FxAvg = FyAvg = TzAvg = PAvg = 0.0;
    }
    if (!std::isfinite(PRigidAvg)) PRigidAvg = 0.0;
    if (!std::isfinite(PDefAvg))   PDefAvg   = 0.0;

    // Froude efficiency
    T thrustPower = std::abs(FxAvg * Vx + FyAvg * Vy) * p.substeps;
    T totalInput  = std::abs(PAvg) * p.substeps;
    T eta = (totalInput > 1e-8) ? thrustPower / std::max(totalInput, T(1e-12)) : 0.0;
    eta = std::min(eta, T(1.0));

    // ---- Body-frame decomposition ----
    // Geometry uses s=0 at the head and s=L at the tail, so forward swimming
    // points toward the head, opposite to the increasing-s spine axis.
    const auto eForward = bodyForwardAxis(theta);
    const auto eLateral = bodyLateralAxis(theta);
    T Uswim    = Vx * eForward[0] + Vy * eForward[1];   // forward swim speed [lu/step]
    T Ulateral = Vx * eLateral[0] + Vy * eLateral[1];  // lateral drift      [lu/step]
    // Net hydrodynamic force projected onto the body forward axis.  Includes
    // drag, pressure, and viscous reactions — this is NOT a separated
    // propulsive thrust component.
    T forwardNetForce = FxAvg * eForward[0] + FyAvg * eForward[1];
    T Flat     = FxAvg * eLateral[0] + FyAvg * eLateral[1];

    // ---- Non-dimensional diagnostics ----
    T absUswim = std::abs(Uswim);
    T bodySpeed = std::sqrt(Vx * Vx + Vy * Vy);
    // Mach number (lattice).
    T machLat = bodySpeed / std::sqrt(T(1) / T(3));
    // Reynolds number: Re = U_swim * total capsule length / nu.
    T Re_inst = absUswim * bodyLength / p.nu();
    // Strouhal number:  St = f * A_pp / |U_swim|
    //   f = eelFreq [Hz],  dtLbm [s/step],  A_pp [lu],  U [lu/step]
    //   St = (f * dtLbm) * A_pp / |U|   (all in lattice units)
    T St_inst = (absUswim > 1e-12)
              ? freqLat * tailPP / absUswim
              : 0.0;
    // Non-dimensional swimming speed:  U* = U_swim / (f * L)
    T Ustar = (absUswim > 1e-12 && freqLat > 1e-12)
            ? absUswim / (freqLat * bodyLength)
            : 0.0;

    // Store history
    const T tFrame = tStart + p.dtAnim;

    history.append({
      tFrame, Vx, Vy, omegaZ,
      FxAvg, FyAvg, TzAvg,
      xCm, yCm, theta,
      PAvg, eta,
      Uswim, Ulateral, forwardNetForce, Flat,
      Re_inst, St_inst, machLat, Ustar,
      meanSlipFrame, maxSlipFrame,
      meanMarkerForceFrame, maxMarkerForceFrame,
      normalizedSlipFrame,
      PRigidAvg, PDefAvg,
      meanResidualSlipFrame, maxResidualSlipFrame,
      meanSoftFluidTorqueFrame, maxSoftFluidTorqueFrame,
      meanSoftAngleStepFrame, maxSoftAngleStepFrame,
      softCouplingResidualFrame, maxSoftCouplingResidualFrame,
      softCouplingItersUsedFrame, maxSoftCouplingItersFrame,
      meanSoftElasticEnergyFrame, maxSoftElasticEnergyFrame,
      meanSoftDampingPowerFrame,
      meanSoftActuatorPowerProxyFrame,
      meanAbsSoftActuatorPowerProxyFrame,
      maxAbsSoftActuatorPowerFrame,
      meanSoftAppliedFluidPowerFrame,
      meanAbsSoftAppliedFluidPowerFrame,
      maxAbsSoftAppliedFluidPowerFrame
    });

    // Console output
    clout << "frame " << frame << "/" << nFrames
          << "  t=" << tFrame
          << "  U*=" << Ustar
          << "  Uswim=" << Uswim
          << "  Ulat=" << Ulateral
          << "  Re=" << Re_inst
          << "  St=" << St_inst
          << "  Ma=" << machLat
          << "  Fnet=" << forwardNetForce
          << "  meanSlip=" << meanSlipFrame
          << "  maxSlip=" << maxSlipFrame
          << "  meanRes=" << meanResidualSlipFrame
          << "  Pdef=" << PDefAvg
          << "  Prigid=" << PRigidAvg
          << "  meanMarkerForce=" << meanMarkerForceFrame
          << "  maxMarkerForce=" << maxMarkerForceFrame
          << "  softTau=" << meanSoftFluidTorqueFrame
          << "  softStep=" << meanSoftAngleStepFrame
          << "  softCoupleRes=" << softCouplingResidualFrame
          << "  softCoupleIter=" << softCouplingItersUsedFrame
          << "  Esoft=" << meanSoftElasticEnergyFrame
          << "  Pdamp=" << meanSoftDampingPowerFrame
          << "  PactProxy=" << meanAbsSoftActuatorPowerProxyFrame
          << "  th=" << (theta * 180.0 / M_PI) << "deg"
          << "  eta=" << eta
          << std::endl;
    if (p.warnMeanSlip > T(0) && meanSlipFrame > p.warnMeanSlip) {
      clout << "WARNING: meanSlip exceeded threshold at frame " << frame
            << "  meanSlip=" << meanSlipFrame
            << "  threshold=" << p.warnMeanSlip << std::endl;
    }
    if (p.warnMaxSlip > T(0) && maxSlipFrame > p.warnMaxSlip) {
      clout << "WARNING: maxSlip exceeded threshold at frame " << frame
            << "  maxSlip=" << maxSlipFrame
            << "  threshold=" << p.warnMaxSlip << std::endl;
    }
    if (p.warnMarkerForce > T(0) && maxMarkerForceFrame > p.warnMarkerForce) {
      clout << "WARNING: marker force exceeded threshold at frame " << frame
            << "  maxMarkerForce=" << maxMarkerForceFrame
            << "  threshold=" << p.warnMarkerForce << std::endl;
    }

    if (softBackboneDynamics &&
        bodyKinematics == BodyKinematics::SoftBackbone &&
        config.softBackboneAbortOnInstability) {
      const T stepLimit = config.softBackboneMaxAngleStep;
      const bool saturatedStep =
        (stepLimit > T(0)) &&
        std::isfinite(maxSoftAngleStepFrame) &&
        maxSoftAngleStepFrame >= T(0.95) * stepLimit;
      softBackboneSaturatedFrameCount =
        saturatedStep ? (softBackboneSaturatedFrameCount + 1) : 0;

      const bool nonFiniteDiagnostics =
        !std::isfinite(meanSlipFrame) ||
        !std::isfinite(maxSlipFrame) ||
        !std::isfinite(meanSoftAngleStepFrame) ||
        !std::isfinite(maxSoftAngleStepFrame) ||
        !std::isfinite(meanSoftFluidTorqueFrame) ||
        !std::isfinite(maxSoftFluidTorqueFrame) ||
        !std::isfinite(softCouplingResidualFrame) ||
        !std::isfinite(maxSoftCouplingResidualFrame) ||
        !std::isfinite(meanSoftElasticEnergyFrame) ||
        !std::isfinite(meanSoftDampingPowerFrame) ||
        !std::isfinite(meanAbsSoftActuatorPowerProxyFrame) ||
        !std::isfinite(meanAbsSoftAppliedFluidPowerFrame);
      const bool excessiveSlip =
        (config.softBackboneAbortMeanSlip > T(0) &&
         meanSlipFrame > config.softBackboneAbortMeanSlip) ||
        (config.softBackboneAbortMaxSlip > T(0) &&
         maxSlipFrame > config.softBackboneAbortMaxSlip);
      const bool repeatedSaturation =
        softBackboneSaturatedFrameCount >=
        config.softBackboneAbortSaturatedFrames;

      if (nonFiniteDiagnostics || excessiveSlip || repeatedSaturation) {
        clout << "ERROR: soft-backbone dynamics instability guard triggered "
              << "at frame " << frame
              << "  t=" << tFrame
              << "  meanSlip=" << meanSlipFrame
              << "  maxSlip=" << maxSlipFrame
              << "  meanSoftFluidTorqueNm=" << meanSoftFluidTorqueFrame
              << "  maxSoftFluidTorqueNm=" << maxSoftFluidTorqueFrame
              << "  meanSoftAngleStep=" << meanSoftAngleStepFrame
              << "  maxSoftAngleStep=" << maxSoftAngleStepFrame
              << "  softCouplingResidual=" << softCouplingResidualFrame
              << "  maxSoftCouplingResidual=" << maxSoftCouplingResidualFrame
              << "  softCouplingItersUsed=" << softCouplingItersUsedFrame
              << "  saturatedFrames="
              << softBackboneSaturatedFrameCount
              << "  limits(meanSlip="
              << config.softBackboneAbortMeanSlip
              << ", maxSlip="
              << config.softBackboneAbortMaxSlip
              << ", saturatedFrames="
              << config.softBackboneAbortSaturatedFrames
              << "). Reduce --softBackboneFluidTorqueScale, increase "
              << "--rampTime, or lower dtLbm via more --substeps."
              << std::endl;
        return 2;
      }
    }

    // Ensure exported fields are synchronized with the latest streamed state.
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Export gating: respect exportInterval and run mode
    const bool doExport = writeSpatialExports && (frame % exportInterval == 0);
    const bool doDiagnostics =
      doExport && exportDiagnostics && (runMode == RunMode::Full);
    const bool doFieldExport =
      doExport && (exportVelocity || exportVorticity || doDiagnostics);
    const bool doBodyExport = doExport && exportBody;

    if (doFieldExport || doBodyExport) {
      LagrangianMarkers frameMarkers;
      if (!buildMarkersChecked(kinematicTimeAt(tFrame), frameMarkers, "export")) {
        return 1;
      }

      // Diagnostic IBM pass: only needed when diagnostics VTI is written.
      // It temporarily rebuilds the instantaneous IBM forcing field for output
      // and is cleared immediately after sampling, so it does not advance or
      // otherwise modify the solver state beyond the intended diagnostic export.
      if (doDiagnostics) {
        ibmStep(p, frameMarkers, xCm, yCm, Vx, Vy, omegaZ, sLattice,
                alphaIBM, ibmIterations, resetEulerForce, ibmRes);
      }

      if (doFieldExport) {
        // Compute body mask once, share across all VTI writers.
        auto bodyMask = computeBodyMask(frameMarkers, nx, ny);
        exportFields.sample(sLattice, bodyMask, doDiagnostics);

        if (exportVelocity) {
          writeVelocityVTI(velDir, frame, vtiWriter, exportFields);
        }
        if (exportVorticity) {
          writeVorticityVTI(vortDir, frame, vtiWriter, exportFields);
        }
        if (doDiagnostics) {
          writeDiagnosticsVTI(diagDir, frame, vtiWriter, exportFields);
        }
      }

      if (doDiagnostics) resetEulerForce();

      if (doBodyExport) {
        writeBodySnapshotCSV(bodyOutDir, frame, frameMarkers);
        writeBodyVTP(bodyOutDir, frame, p, frameMarkers);
      }
    }
  }

  // ============================================================
  //  Save ParaView .pvd time-series files (one per field group)
  // ============================================================
  if (writeSpatialExports) {
    bool wrotePvd = false;
    if (exportVelocity) {
      std::string pvdPath = writePVD(runOutDir, "eel3dof_velocity.pvd", "vtkData/velocity/", "eel3dof_velocity_", ".vti", nFrames, exportInterval, p.dtAnim);
      clout << "Saved PVD: " << pvdPath << std::endl;
      wrotePvd = true;
    }
    if (exportVorticity) {
      std::string pvdPath = writePVD(runOutDir, "eel3dof_vorticity.pvd", "vtkData/vorticity/", "eel3dof_vorticity_", ".vti", nFrames, exportInterval, p.dtAnim);
      clout << "Saved PVD: " << pvdPath << std::endl;
      wrotePvd = true;
    }
    if (runMode == RunMode::Full && exportDiagnostics) {
      std::string pvdPath = writePVD(runOutDir, "eel3dof_diagnostics.pvd", "vtkData/diagnostics/", "eel3dof_diagnostics_", ".vti", nFrames, exportInterval, p.dtAnim);
      clout << "Saved PVD: " << pvdPath << std::endl;
      wrotePvd = true;
    }
    if (exportBody) {
      std::string pvdPath = writePVD(runOutDir, "eel3dof_body.pvd", "bodyData/", "eel3dof_body_", ".vtp", nFrames, exportInterval, p.dtAnim);
      clout << "Saved PVD: " << pvdPath << std::endl;
      wrotePvd = true;
    }
    if (!wrotePvd) {
      clout << "Spatial exports disabled by export flags." << std::endl;
    }
  } else {
    clout << "Verification mode: skipped VTK/body snapshot exports." << std::endl;
  }

  // ============================================================
  //  Append one-run aspect-ratio summary
  // ============================================================
  SteadySummary steady = computeSteadySummary(
    history.histT, history.histUswim, history.histPower,
    history.histForwardNetForce, history.histFlat,
    history.histRe, history.histSt, history.histUstar, history.histEta,
    history.histMeanSlip, history.histMaxSlip,
    history.histMeanMarkerForce, history.histMaxMarkerForce,
    history.histNormalizedSlip,
    history.histPowerRigid, history.histPowerDef,
    history.histMeanResidualSlip, history.histMaxResidualSlip,
    tCut, mass);

  if (steady.nSamples == 0) {
    clout << "WARNING: no steady-state samples found after tCut=" << tCut
          << ". Increase Ttotal or lower --tCut for summary metrics." << std::endl;
  }

  auto meanAfterTCut = [&](const std::vector<T>& values) -> T {
    T sum = T(0);
    int n = 0;
    const size_t nLimit = std::min(history.histT.size(), values.size());
    for (size_t i = 0; i < nLimit; ++i) {
      if (!(history.histT[i] > tCut) || !std::isfinite(values[i])) continue;
      sum += values[i];
      ++n;
    }
    return (n > 0) ? (sum / T(n)) : T(0);
  };
  auto maxAfterTCut = [&](const std::vector<T>& values) -> T {
    T maxValue = T(0);
    const size_t nLimit = std::min(history.histT.size(), values.size());
    for (size_t i = 0; i < nLimit; ++i) {
      if (!(history.histT[i] > tCut) || !std::isfinite(values[i])) continue;
      maxValue = std::max(maxValue, std::abs(values[i]));
    }
    return maxValue;
  };
  const T steadyMeanSoftFluidTorqueNm =
    meanAfterTCut(history.histMeanSoftFluidTorqueNm);
  const T steadyMaxSoftFluidTorqueNm =
    maxAfterTCut(history.histMaxSoftFluidTorqueNm);
  const T steadyMeanSoftAngleStep = meanAfterTCut(history.histMeanSoftAngleStep);
  const T steadyMaxSoftAngleStep = maxAfterTCut(history.histMaxSoftAngleStep);
  const T steadyMeanSoftCouplingResidual =
    meanAfterTCut(history.histSoftCouplingResidual);
  const T steadyMaxSoftCouplingResidual =
    maxAfterTCut(history.histMaxSoftCouplingResidual);
  const T steadyMeanSoftCouplingItersUsed =
    meanAfterTCut(history.histSoftCouplingItersUsed);
  const T steadyMaxSoftCouplingItersUsed =
    maxAfterTCut(history.histMaxSoftCouplingItersUsed);
  const T steadyMeanSoftElasticEnergyJ =
    meanAfterTCut(history.histSoftElasticEnergyJ);
  const T steadyMaxSoftElasticEnergyJ =
    maxAfterTCut(history.histMaxSoftElasticEnergyJ);
  const T steadyMeanSoftDampingPowerW =
    meanAfterTCut(history.histSoftDampingPowerW);
  const T steadyMeanSoftActuatorPowerProxyW =
    meanAfterTCut(history.histSoftActuatorPowerProxyW);
  const T steadyMeanAbsSoftActuatorPowerProxyW =
    meanAfterTCut(history.histAbsSoftActuatorPowerProxyW);
  const T steadyMaxAbsSoftActuatorPowerProxyW =
    maxAfterTCut(history.histMaxAbsSoftActuatorPowerProxyW);
  const T steadyMeanSoftAppliedFluidPowerW =
    meanAfterTCut(history.histSoftAppliedFluidPowerW);
  const T steadyMeanAbsSoftAppliedFluidPowerW =
    meanAfterTCut(history.histAbsSoftAppliedFluidPowerW);
  const T steadyMaxAbsSoftAppliedFluidPowerW =
    maxAfterTCut(history.histMaxAbsSoftAppliedFluidPowerW);

  const T dxM = (p.totalGeometricLengthLU() > T(0))
              ? p.physicalBodyLengthM / p.totalGeometricLengthLU()
              : T(0);
  const T steadyMeanUPhysicalMps =
    (dxM > T(0) && dtLbm > T(0)) ? steady.meanU * dxM / dtLbm : T(0);
  const T softBodyMassKg =
    (rodSection.valid && rodSection.massPerLengthKgM > T(0))
      ? rodSection.massPerLengthKgM * rodSection.physicalLengthM
      : T(0);
  const T cotSoftActuatorProxySI =
    (softBodyMassKg > T(1e-12) && steadyMeanUPhysicalMps > T(1e-12))
      ? steadyMeanAbsSoftActuatorPowerProxyW
        / (softBodyMassKg * steadyMeanUPhysicalMps)
      : T(0);

  const int recommendedSteadySamples = recommendedSteadySampleCount(p);
  const bool enoughSteadySamples = (steady.nSamples >= recommendedSteadySamples);

  // ----------------------------------------------------------------
  //  Cycle-averaged steady-state diagnostics
  // ----------------------------------------------------------------
  //  period = 1 / eelFreq (physical seconds, same units as history time
  //  samples).  Bins start at the first sample with t > tCut.  Convergence is
  //  reported over the last min(5, nCycles) complete cycles.
  const T cyclePeriod = (p.eelFreq > T(1e-12)) ? (T(1) / p.eelFreq) : T(0);
  std::vector<CycleAverage> cycles = computeCycleAverages(
    history.histT, history.histUswim, history.histPower,
    history.histForwardNetForce,
    history.histRe, history.histSt, history.histUstar,
    history.histMeanSlip, history.histNormalizedSlip,
    tCut, cyclePeriod, mass);
  const CycleConvergence cycleConv = computeCycleConvergence(cycles);

  SummaryCsvInputs summaryInput;
  summaryInput.p = p;
  summaryInput.steady = steady;
  summaryInput.cycleConv = cycleConv;
  summaryInput.simCase = simCase;
  summaryInput.warmupMode = warmupMode;
  summaryInput.wallBoundary = wallBoundary;
  summaryInput.bodyKinematics = bodyKinematics;
  summaryInput.initialPlacementClamped = clampAudit.initialPlacementClamped;
  summaryInput.initialPlacementClampCount = clampAudit.initialPlacementClampCount;
  summaryInput.runtimeDomainClampHit = clampAudit.runtimeDomainClampHit;
  summaryInput.runtimeDomainClampCount = clampAudit.runtimeDomainClampCount;
  summaryInput.nx = nx;
  summaryInput.ny = ny;
  summaryInput.tCut = tCut;
  summaryInput.bodyCenterlineLength = bodyCenterlineLength;
  summaryInput.bodyLength = bodyLength;
  summaryInput.bodyWidth = bodyWidth;
  summaryInput.mass = mass;
  summaryInput.Ibody = Ibody;
  summaryInput.alphaIBM = alphaIBM;
  summaryInput.ibmIterations = ibmIterations;
  summaryInput.softBackboneDynamics = softBackboneDynamics;
  summaryInput.softBackboneRelaxationTime = config.softBackboneRelaxationTime;
  summaryInput.softBackboneFluidTorqueScale =
    config.softBackboneFluidTorqueScale;
  summaryInput.softBackboneFluidTorqueFilterTime =
    config.softBackboneFluidTorqueFilterTime;
  summaryInput.softBackboneAddedMassFrac =
    config.softBackboneAddedMassFrac;
  summaryInput.softBackboneMaxAngleStep = config.softBackboneMaxAngleStep;
  summaryInput.softBackboneCouplingIterations =
    config.softBackboneCouplingIterations;
  summaryInput.softBackboneCouplingRelaxation =
    config.softBackboneCouplingRelaxation;
  summaryInput.softBackboneCouplingTolerance =
    config.softBackboneCouplingTolerance;
  summaryInput.softBackboneLoadProjection =
    config.softBackboneLoadProjection;
  summaryInput.meanSoftFluidTorqueNm = steadyMeanSoftFluidTorqueNm;
  summaryInput.maxSoftFluidTorqueNm = steadyMaxSoftFluidTorqueNm;
  summaryInput.meanSoftAngleStep = steadyMeanSoftAngleStep;
  summaryInput.maxSoftAngleStep = steadyMaxSoftAngleStep;
  summaryInput.meanSoftCouplingResidual =
    steadyMeanSoftCouplingResidual;
  summaryInput.maxSoftCouplingResidual =
    steadyMaxSoftCouplingResidual;
  summaryInput.meanSoftCouplingItersUsed =
    steadyMeanSoftCouplingItersUsed;
  summaryInput.maxSoftCouplingItersUsed =
    steadyMaxSoftCouplingItersUsed;
  summaryInput.meanSoftElasticEnergyJ =
    steadyMeanSoftElasticEnergyJ;
  summaryInput.maxSoftElasticEnergyJ =
    steadyMaxSoftElasticEnergyJ;
  summaryInput.meanSoftDampingPowerW =
    steadyMeanSoftDampingPowerW;
  summaryInput.meanSoftActuatorPowerProxyW =
    steadyMeanSoftActuatorPowerProxyW;
  summaryInput.meanAbsSoftActuatorPowerProxyW =
    steadyMeanAbsSoftActuatorPowerProxyW;
  summaryInput.maxAbsSoftActuatorPowerProxyW =
    steadyMaxAbsSoftActuatorPowerProxyW;
  summaryInput.meanSoftAppliedFluidPowerW =
    steadyMeanSoftAppliedFluidPowerW;
  summaryInput.meanAbsSoftAppliedFluidPowerW =
    steadyMeanAbsSoftAppliedFluidPowerW;
  summaryInput.maxAbsSoftAppliedFluidPowerW =
    steadyMaxAbsSoftAppliedFluidPowerW;
  summaryInput.softBodyMassKg = softBodyMassKg;
  summaryInput.meanUPhysicalMps = steadyMeanUPhysicalMps;
  summaryInput.cotSoftActuatorProxySI = cotSoftActuatorProxySI;

  const CsvWriteResult arCsvResult = appendArSummaryCsv(summaryCsv, summaryInput);
  if (arCsvResult.resolution.headerMismatch) {
    clout << "WARNING: existing CSV header in " << summaryCsv
          << " does not match the current schema. Writing this run to "
          << arCsvResult.resolution.fallbackPath << " instead." << std::endl;
  }
  if (arCsvResult.opened) {
    clout << "Saved AR summary row: " << arCsvResult.resolution.path << std::endl;
  } else {
    clout << "WARNING: could not open summary CSV: " << arCsvResult.resolution.path << std::endl;
  }

  VerificationCsvInputs verificationInput;
  verificationInput.p = p;
  verificationInput.steady = steady;
  verificationInput.cycleConv = cycleConv;
  verificationInput.runId = runId;
  verificationInput.simCase = simCase;
  verificationInput.studyMode = studyMode;
  verificationInput.warmupMode = warmupMode;
  verificationInput.wallBoundary = wallBoundary;
  verificationInput.bodyKinematics = bodyKinematics;
  verificationInput.initialPlacementClamped = clampAudit.initialPlacementClamped;
  verificationInput.initialPlacementClampCount = clampAudit.initialPlacementClampCount;
  verificationInput.runtimeDomainClampHit = clampAudit.runtimeDomainClampHit;
  verificationInput.runtimeDomainClampCount = clampAudit.runtimeDomainClampCount;
  verificationInput.nx = nx;
  verificationInput.ny = ny;
  verificationInput.dtLbm = dtLbm;
  verificationInput.tCut = tCut;
  verificationInput.bodyCenterlineLength = bodyCenterlineLength;
  verificationInput.bodyLength = bodyLength;
  verificationInput.bodyWidth = bodyWidth;
  verificationInput.mass = mass;
  verificationInput.Ibody = Ibody;
  verificationInput.alphaIBM = alphaIBM;
  verificationInput.ibmIterations = ibmIterations;
  verificationInput.softBackboneDynamics = softBackboneDynamics;
  verificationInput.softBackboneRelaxationTime =
    config.softBackboneRelaxationTime;
  verificationInput.softBackboneFluidTorqueScale =
    config.softBackboneFluidTorqueScale;
  verificationInput.softBackboneFluidTorqueFilterTime =
    config.softBackboneFluidTorqueFilterTime;
  verificationInput.softBackboneAddedMassFrac =
    config.softBackboneAddedMassFrac;
  verificationInput.softBackboneMaxAngleStep = config.softBackboneMaxAngleStep;
  verificationInput.softBackboneCouplingIterations =
    config.softBackboneCouplingIterations;
  verificationInput.softBackboneCouplingRelaxation =
    config.softBackboneCouplingRelaxation;
  verificationInput.softBackboneCouplingTolerance =
    config.softBackboneCouplingTolerance;
  verificationInput.softBackboneLoadProjection =
    config.softBackboneLoadProjection;
  verificationInput.meanSoftFluidTorqueNm = steadyMeanSoftFluidTorqueNm;
  verificationInput.maxSoftFluidTorqueNm = steadyMaxSoftFluidTorqueNm;
  verificationInput.meanSoftAngleStep = steadyMeanSoftAngleStep;
  verificationInput.maxSoftAngleStep = steadyMaxSoftAngleStep;
  verificationInput.meanSoftCouplingResidual =
    steadyMeanSoftCouplingResidual;
  verificationInput.maxSoftCouplingResidual =
    steadyMaxSoftCouplingResidual;
  verificationInput.meanSoftCouplingItersUsed =
    steadyMeanSoftCouplingItersUsed;
  verificationInput.maxSoftCouplingItersUsed =
    steadyMaxSoftCouplingItersUsed;
  verificationInput.meanSoftElasticEnergyJ =
    steadyMeanSoftElasticEnergyJ;
  verificationInput.maxSoftElasticEnergyJ =
    steadyMaxSoftElasticEnergyJ;
  verificationInput.meanSoftDampingPowerW =
    steadyMeanSoftDampingPowerW;
  verificationInput.meanSoftActuatorPowerProxyW =
    steadyMeanSoftActuatorPowerProxyW;
  verificationInput.meanAbsSoftActuatorPowerProxyW =
    steadyMeanAbsSoftActuatorPowerProxyW;
  verificationInput.maxAbsSoftActuatorPowerProxyW =
    steadyMaxAbsSoftActuatorPowerProxyW;
  verificationInput.meanSoftAppliedFluidPowerW =
    steadyMeanSoftAppliedFluidPowerW;
  verificationInput.meanAbsSoftAppliedFluidPowerW =
    steadyMeanAbsSoftAppliedFluidPowerW;
  verificationInput.maxAbsSoftAppliedFluidPowerW =
    steadyMaxAbsSoftAppliedFluidPowerW;
  verificationInput.softBodyMassKg = softBodyMassKg;
  verificationInput.meanUPhysicalMps = steadyMeanUPhysicalMps;
  verificationInput.cotSoftActuatorProxySI = cotSoftActuatorProxySI;

  const CsvWriteResult sensitivityCsvResult = appendVerificationSummaryCsv(sensitivityCsv, verificationInput);
  if (sensitivityCsvResult.resolution.headerMismatch) {
    clout << "WARNING: existing CSV header in " << sensitivityCsv
          << " does not match the current schema. Writing this run to "
          << sensitivityCsvResult.resolution.fallbackPath << " instead." << std::endl;
  }
  if (sensitivityCsvResult.opened) {
    clout << "Saved verification summary row: " << sensitivityCsvResult.resolution.path << std::endl;
  } else {
    clout << "WARNING: could not open sensitivity CSV: " << sensitivityCsvResult.resolution.path << std::endl;
  }

  clout << "==== Verification Summary ====" << std::endl;
  clout << "run_id=" << runId
        << "  simCase=" << simulationCaseName(simCase)
        << "  grid=" << nx << "x" << ny
        << "  tau=" << p.tau
        << "  substeps=" << p.substeps << std::endl;
  clout << "meanU=" << steady.meanU
        << "  meanP=" << steady.meanP
        << "  CoT=" << steady.CoT
        << "  hydroCost=" << steady.hydroCost
        << "  transportEfficiencyProxy=" << steady.transportEfficiencyProxy
        << "  meanForwardNetForce=" << steady.meanForwardNetForce
        << "  meanAbsForwardNetForce=" << steady.meanAbsForwardNetForce
        << "  meanSignedForwardNetForce=" << steady.meanSignedForwardNetForce
        << "  etaNetForceDiagnostic=" << steady.etaNetForceDiagnostic
        << "  meanSlip=" << steady.meanSlip
        << "  meanNormalizedSlip=" << steady.meanNormalizedSlip
        << "  meanMarkerForce=" << steady.meanMarkerForce << std::endl;
  clout << "NOTE: meanForwardNetForce is the mean net hydrodynamic force "
        << "projected onto the forward axis; for self-propelled steady "
        << "swimming it can approach zero or change sign and is NOT a "
        << "separated propulsive thrust.  etaNetForceDiagnostic is therefore "
        << "diagnostic-only — rank cases by CoT, hydroCost, "
        << "transportEfficiencyProxy, and meanU." << std::endl;
  clout << "cyclePeriod=" << cyclePeriod
        << "  nSteadyCycles=" << cycleConv.nSteadyCycles
        << "  cycleMeanUswim=" << cycleConv.cycleMeanUswim
        << "  cycleCvUswim=" << cycleConv.cycleCvUswim
        << "  cycleMeanUstar=" << cycleConv.cycleMeanUstar
        << "  cycleCvUstar=" << cycleConv.cycleCvUstar
        << "  cycleMeanCoT=" << cycleConv.cycleMeanCoT
        << "  cycleCvCoT=" << cycleConv.cycleCvCoT
        << "  cycleMeanHydroCost=" << cycleConv.cycleMeanHydroCost
        << "  cycleCvHydroCost=" << cycleConv.cycleCvHydroCost
        << "  cycleMeanNormalizedSlip=" << cycleConv.cycleMeanNormalizedSlip
        << "  cycleCvNormalizedSlip=" << cycleConv.cycleCvNormalizedSlip
        << "  cycleConverged=" << (cycleConv.cycleConverged ? "yes" : "no")
        << std::endl;
  clout << "runtimeDomainClampHit=" << (clampAudit.runtimeDomainClampHit ? 1 : 0)
        << "  runtimeDomainClampCount=" << clampAudit.runtimeDomainClampCount
        << std::endl;
  clout << "initialPlacementClamped=" << (clampAudit.initialPlacementClamped ? 1 : 0)
        << "  initialPlacementClampCount=" << clampAudit.initialPlacementClampCount
        << "  (one-shot setup nudge from initialPositionFactor / rightMargin; "
        << "not a runtime issue)" << std::endl;
  clout << "warmupMode=" << warmupModeName(warmupMode)
        << "  waveDirection=" << waveDirectionName(p.waveDirection)
        << "  bodyKinematics=" << bodyKinematicsName(bodyKinematics)
        << "  eelFreq=" << p.eelFreq
        << "  eelA0=" << p.eelA0 << std::endl;
  clout << "IBM(v7): alphaIBM=" << alphaIBM
        << "  ibmIterations=" << ibmIterations
        << "  meanResidualSlip=" << steady.meanResidualSlip
        << "  maxResidualSlip=" << steady.maxResidualSlip << std::endl;
  if (bodyKinematics == BodyKinematics::SoftBackbone) {
    clout << "Soft backbone dynamics: enabled="
          << (softBackboneDynamics ? 1 : 0)
          << "  fluidTorqueScale=" << config.softBackboneFluidTorqueScale
          << "  fluidTorqueFilterTime="
          << config.softBackboneFluidTorqueFilterTime
          << "  maxAngleStepLimit=" << config.softBackboneMaxAngleStep
          << "  couplingIterations="
          << config.softBackboneCouplingIterations
          << "  couplingRelaxation="
          << config.softBackboneCouplingRelaxation
          << "  couplingTolerance="
          << config.softBackboneCouplingTolerance
          << "  loadProjection="
          << softBackboneLoadProjectionName(
               config.softBackboneLoadProjection)
          << "  meanSoftFluidTorqueNm=" << steadyMeanSoftFluidTorqueNm
          << "  maxSoftFluidTorqueNm=" << steadyMaxSoftFluidTorqueNm
          << "  meanSoftAngleStep=" << steadyMeanSoftAngleStep
          << "  maxSoftAngleStep=" << steadyMaxSoftAngleStep
          << "  meanSoftCouplingResidual="
          << steadyMeanSoftCouplingResidual
          << "  maxSoftCouplingResidual="
          << steadyMaxSoftCouplingResidual
          << "  meanSoftCouplingItersUsed="
          << steadyMeanSoftCouplingItersUsed
          << "  maxSoftCouplingItersUsed="
          << steadyMaxSoftCouplingItersUsed
          << "  meanSoftElasticEnergyJ="
          << steadyMeanSoftElasticEnergyJ
          << "  maxSoftElasticEnergyJ="
          << steadyMaxSoftElasticEnergyJ
          << "  meanSoftDampingPowerW="
          << steadyMeanSoftDampingPowerW
          << "  meanAbsSoftActuatorPowerProxyW="
          << steadyMeanAbsSoftActuatorPowerProxyW
          << "  meanAbsSoftAppliedFluidPowerW="
          << steadyMeanAbsSoftAppliedFluidPowerW
          << "  meanUPhysicalMps="
          << steadyMeanUPhysicalMps
          << "  cotSoftActuatorProxySI="
          << cotSoftActuatorProxySI << std::endl;
  }
  clout << "Power decomposition (steady window): "
        << "meanConstraintPower=" << steady.meanConstraintPower
        << "  meanRigidBodyPower=" << steady.meanRigidBodyPower
        << "  meanDeformationPower=" << steady.meanDeformationPower
        << "  meanAbsConstraintPower=" << steady.meanAbsConstraintPower
        << "  meanAbsRigidBodyPower=" << steady.meanAbsRigidBodyPower
        << "  meanAbsDeformationPower=" << steady.meanAbsDeformationPower
        << "  transportEfficiencyDef=" << steady.transportEfficiencyDef
        << "  cotDef=" << steady.cotDef << std::endl;
  {
    const T sumDecomp = steady.meanRigidBodyPower + steady.meanDeformationPower;
    const T denom = std::max(std::abs(steady.meanConstraintPower), T(1e-30));
    const T relErr = std::abs(sumDecomp - steady.meanConstraintPower) / denom;
    clout << "Power identity check: meanConstraintPower vs "
          << "meanRigidBodyPower+meanDeformationPower  ("
          << steady.meanConstraintPower << " vs " << sumDecomp
          << ", rel.err=" << relErr << ")" << std::endl;
  }
  clout << "steadySamples=" << steady.nSamples
        << "  recommended>=" << recommendedSteadySamples
        << "  enough=" << (enoughSteadySamples ? "yes" : "no")
        << std::endl;

  // ============================================================
  //  Save time-series data
  // ============================================================
  {
    const std::string fname = writeHistoryCsv(runOutDir, history.csvData());
    clout << "Saved history: " << fname << std::endl;
  }

  // ---- Write marker positions for last frame (body outline) ----
  {
    const std::string fname = writeFinalBodyCsv(runOutDir, markers);
    clout << "Saved final body: " << fname << std::endl;
  }

  clout << "Simulation complete." << std::endl;
  return 0;
}
