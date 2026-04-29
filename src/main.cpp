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
#include "physics/diagnostics.hpp"
#include "physics/rigid_body.hpp"
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

  RunConfig config = parseCommandLine(argc, argv, baseLogDir);
  EelParams& p = config.p;
  RunMode& runMode = config.runMode;
  StudyMode& studyMode = config.studyMode;
  SimulationCase& simCase = config.simCase;
  WarmupMode& warmupMode = config.warmupMode;
  GaitNormalization& gaitNormalization = config.gaitNormalization;
  T& tailAmpRatioTarget = config.tailAmpRatioTarget;
  T& targetSt = config.targetSt;
  T& referenceU = config.referenceU;
  T& alphaIBM = config.alphaIBM;
  int& ibmIterations = config.ibmIterations;
  bool& legacyKappaInputUsed = config.legacyKappaInputUsed;
  T& legacyKappaInputValue = config.legacyKappaInputValue;
  bool& rawGeometryOverride = config.rawGeometryOverride;
  bool& aspectGeometryOverride = config.aspectGeometryOverride;
  T& tCut = config.tCut;
  std::string& summaryCsv = config.summaryCsv;
  std::string& sensitivityCsv = config.sensitivityCsv;
  std::string& runTag = config.runTag;

  if (legacyKappaInputUsed) {
    clout << "Deprecated: --kappa mapped to --alphaIBM (kappa="
          << legacyKappaInputValue << " -> alphaIBM=" << alphaIBM
          << "). Prefer --alphaIBM=<value> in new scripts." << std::endl;
  }

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

  // ----------------------------------------------------------------
  //  Gait normalization for AR sweeps
  // ----------------------------------------------------------------
  //  Applied AFTER geometry has been resolved.  Preserves restTime /
  //  rampTime / lambda / substeps unchanged; only eelA0 (tailAmpRatio)
  //  or eelFreq (targetSt) is rescaled.  The pre-normalization values are
  //  recorded so they can be reported alongside the effective values.
  const T preNormEelFreq = p.eelFreq;
  const T preNormEelA0   = p.eelA0;
  GaitNormalization effectiveGaitNormalization = gaitNormalization;
  T resolvedTailAmpRatioTarget = tailAmpRatioTarget;
  if (gaitNormalization == GaitNormalization::TailAmpRatio) {
    // Baseline = default EelParams geometry & gait.  This gives a single
    // well-defined target ratio that AR sweeps can be normalized against
    // even when the user does not pass --tailAmpRatioTarget explicitly.
    EelParams pBaseline;
    if (!updateGeometryFromAspectRatio(pBaseline)) {
      clout << "WARNING: tailAmpRatio normalization could not derive a "
            << "default baseline; falling back to gaitNormalization=fixed."
            << std::endl;
      effectiveGaitNormalization = GaitNormalization::Fixed;
    } else {
      const T baselineRatio =
        pBaseline.tailAmplitudeLU() / pBaseline.totalGeometricLengthLU();
      const T targetRatio = (tailAmpRatioTarget > T(0))
                          ? tailAmpRatioTarget : baselineRatio;
      resolvedTailAmpRatioTarget = targetRatio;
      // tailAmplitudeLU = eelA0 * 2.3 * eelScale,
      // totalLength    = eelScale + 2*bodyRadius
      // => eelA0 = targetRatio * totalLength / (2.3 * eelScale)
      const T totalLength = p.totalGeometricLengthLU();
      const T denom = T(2.3) * p.eelScale;
      if (denom > T(1e-12)) {
        p.eelA0 = targetRatio * totalLength / denom;
      } else {
        clout << "WARNING: tailAmpRatio normalization aborted (eelScale=0); "
              << "falling back to gaitNormalization=fixed." << std::endl;
        effectiveGaitNormalization = GaitNormalization::Fixed;
      }
    }
  } else if (gaitNormalization == GaitNormalization::TargetSt) {
    // Strouhal in lattice units:
    //   St = (eelFreq[Hz] * dtLbm[s/step]) * tailPP[lu] / U[lu/step]
    // referenceU is given in lattice units (lu/step) to match the reported
    // Uswim.  Solve for eelFreq[Hz]:
    //   eelFreq = targetSt * U / (dtLbm * tailPP)
    const T tailPP = p.tailPeakToPeakLU();
    const T dtLbm  = p.dtLbm();
    if (!(targetSt > T(0)) || !(referenceU > T(0)) ||
        !(tailPP > T(1e-12)) || !(dtLbm > T(1e-12))) {
      clout << "WARNING: targetSt normalization needs --targetSt>0 and "
            << "--referenceU>0 with a finite tailPP and dtLbm; falling back "
            << "to gaitNormalization=fixed." << std::endl;
      effectiveGaitNormalization = GaitNormalization::Fixed;
    } else {
      p.eelFreq = targetSt * referenceU / (dtLbm * tailPP);
    }
  }
  const T effectiveEelFreq = p.eelFreq;
  const T effectiveEelA0   = p.eelA0;

  if (tCut < T(0)) {
    const T periodsForAverage = (p.eelFreq > T(1e-12)) ? T(5) / p.eelFreq : T(3);
    tCut = p.restTime + p.rampTime + std::max(T(3), periodsForAverage);
  }

  const std::string runsRootDir = baseLogDir + "runs/";
  if (!eelEnsureDirectoryTree(runsRootDir)) {
    clout << "ERROR: could not create run output root: " << runsRootDir << std::endl;
    return 1;
  }
  const std::string baseRunId = makeBaseRunId(p, simCase, runTag);
  const std::string runId = makeUniqueRunId(baseRunId, runsRootDir);
  const std::string runOutDir = runsRootDir + runId + "/";
  const std::string runVtkDir = runOutDir + "vtkData/";
  const std::string bodyOutDir = runOutDir + "bodyData/";
  const std::string velDir  = runVtkDir + "velocity/";
  const std::string vortDir = runVtkDir + "vorticity/";
  const std::string diagDir = runVtkDir + "diagnostics/";
  if (!eelEnsureDirectoryTree(bodyOutDir) ||
      !eelEnsureDirectoryTree(velDir) ||
      !eelEnsureDirectoryTree(vortDir) ||
      !eelEnsureDirectoryTree(diagDir)) {
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
  const bool fixedInflow = (simCase == SimulationCase::FixedInflow);

  clout << "==== Self-propelled eel 3-DOF (OpenLB-native fluid + custom IBM) ====" << std::endl;
  clout << "Run ID: " << runId << std::endl;
  clout << "Run output: " << runOutDir << std::endl;
  clout << "Study mode: " << studyModeName(studyMode) << std::endl;
  clout << "Sensitivity summary: " << sensitivityCsv << std::endl;
  clout << "Grid: " << nx << " x " << ny << "  tau=" << p.tau
        << "  nu=" << p.nu() << std::endl;
  clout << "Case: " << simulationCaseName(simCase)
        << "  inflowU=" << (fixedInflow ? p.inflowVelocity : T(0)) << std::endl;
  clout << "AR study: useAspectRatioGeometry=" << (p.useAspectRatioGeometry ? "true" : "false")
        << "  AR=" << p.aspectRatio
        << "  bodyAreaTarget=" << p.bodyAreaTarget
        << "  tCut=" << tCut << std::endl;
  if (simCase != SimulationCase::SurgeOnly) {
    clout << "NOTE: surge_only is the preferred first optimization mode for AR studies; "
          << "other cases are retained for validation/debugging." << std::endl;
  } else {
    clout << "Mode note: surge_only advances only the streamwise rigid-body DOF; "
          << "heave/yaw are excluded by model definition." << std::endl;
  }
  clout << "Eel: centerlineL=" << p.centerlineLengthLU()
        << " lu  totalL=" << p.totalGeometricLengthLU()
        << " lu  W=" << p.bodyWidthLU()
        << " lu  freq=" << p.eelFreq << "  lambda=" << p.eelLambda << std::endl;
  clout << "Warmup: mode=" << warmupModeName(warmupMode)
        << "  nWarmup=" << p.nWarmup << std::endl;
  clout << "Gait normalization: requested=" << gaitNormalizationName(gaitNormalization)
        << "  effective=" << gaitNormalizationName(effectiveGaitNormalization)
        << "  preNorm(eelFreq=" << preNormEelFreq
        << ", eelA0=" << preNormEelA0 << ")"
        << "  effective(eelFreq=" << effectiveEelFreq
        << ", eelA0=" << effectiveEelA0 << ")";
  if (effectiveGaitNormalization == GaitNormalization::TailAmpRatio) {
    clout << "  tailAmpRatioTarget=" << resolvedTailAmpRatioTarget;
  } else if (effectiveGaitNormalization == GaitNormalization::TargetSt) {
    clout << "  targetSt=" << targetSt
          << "  referenceU=" << referenceU << " (lu/step)";
  }
  clout << std::endl;
  clout << "Time: Ttotal=" << p.Ttotal << "  dt_anim=" << p.dtAnim
        << "  substeps=" << p.substeps << "  dt_lbm=" << dtLbm << std::endl;
  clout << "OpenMP threads: " << omp_get_max_threads()
        << "  mode=" << (runMode == RunMode::Preview ? "preview" :
                         runMode == RunMode::Standard ? "standard" : "full")
        << "  exportInterval=" << exportInterval << std::endl;
  if (verificationMode) {
    clout << "Verification mode keeps the same solver physics/update order and "
          << "suppresses VTK/body exports to make parameter sweeps lighter."
          << std::endl;
  }
  clout << "IBM: alphaIBM=" << alphaIBM
        << "  ibmIterations=" << ibmIterations
        << "  legacyKappaInput=" << (legacyKappaInputUsed ? 1 : 0)
        << " (forcing law: F = alphaIBM * rho_local * (Ud - U_interp) / dt;"
        << " rho_local=1.0, dt=1.0 in lattice units)"
        << std::endl;
  clout << "IBM legacy params: kappa(stored)=" << p.kappa
        << "  nIbmIters(requested)=" << p.nIbmIters
        << "  spongeWidth=" << p.spongeWidth
        << "  spongeStrength=" << p.spongeStrength << std::endl;
  clout << "IBM warnings: meanSlip>" << p.warnMeanSlip
        << "  maxSlip>" << p.warnMaxSlip
        << "  maxMarkerForce>" << p.warnMarkerForce << std::endl;
  if (p.nIbmIters != 1) {
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

  // Channel layout follows the Python vortex-street validation case:
  // top/bottom are no-slip walls; fixed_inflow uses a left velocity inlet
  // and right pressure outlet, while swimming cases keep both sides open.
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
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  if (fixedInflow) {
    boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  } else {
    boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 3);
  }
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
  sLattice.template setParameter<descriptors::OMEGA>(p.omega_lbm());

  AnalyticalConst2D<T,T> rho0(1.0);
  AnalyticalConst2D<T,T> u0(fixedInflow ? p.inflowVelocity : T(0), 0.0);
  AnalyticalConst2D<T,T> f0(0.0, 0.0);
  auto fluidAndBoundary = superGeometry.getMaterialIndicator({1, 2, 3, 4});
  sLattice.iniEquilibrium(fluidAndBoundary, rho0, u0);
  sLattice.defineRhoU(fluidAndBoundary, rho0, u0);
  if (fixedInflow) {
    AnalyticalConst2D<T,T> uIn(p.inflowVelocity, 0.0);
    sLattice.defineU(superGeometry, 3, uIn);
  }
  fields::set<descriptors::FORCE>(sLattice, fluidAndBoundary, f0);
  sLattice.initialize();
  if (fixedInflow) {
    sLattice.template setProcessingContext<
      Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
  }

  // ---- Sponge zone near open (left/right) boundaries ----
  //
  // The previous solver blended raw populations toward f_eq.  In the OpenLB
  // workflow we keep this as a lattice FORCE-field damping layer, so the
  // damping is applied by ForcedBGKdynamics instead of direct population edits.
  auto resetEulerForce = [&]() {
    resetEulerForceField(p, sLattice, fixedInflow);
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

  // Body inertial properties (neutrally buoyant).  The centerline length is
  // used by the gait construction; the capsule tip-to-tip length is the
  // characteristic body length for mass/inertia/Re diagnostics.
  BodyInertia bodyInertia = computeBodyInertia(p, 1.0);
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

  // Validation-friendly case semantics:
  //  fixed_inflow      : body held fixed and undulation disabled, so inlet/outlet
  //                      and wall behavior can be checked against a stationary body.
  //  fixed_undulation  : body pose held fixed while the prescribed undulation remains
  //                      active, isolating IBM enforcement and wake generation.
  //  surge_only        : true reduced-order self-propelled model; only streamwise
  //                      translation is a dynamic state. Heave/yaw are not evolved.
  auto deformationAmpScale = [&]() -> T {
    return (simCase == SimulationCase::FixedInflow) ? T(0) : T(1);
  };

  const bool fixedPoseMode =
    (simCase == SimulationCase::FixedInflow ||
     simCase == SimulationCase::FixedUndulation);
  const bool surgeOnlyMode = (simCase == SimulationCase::SurgeOnly);
  const bool full3DofMode  = (simCase == SimulationCase::Full3Dof);

  auto applyModeDefinitionState = [&]() {
    enforceModeDefinitionState(simCase, {xCmRef, yCmRef, thetaRef}, bodyState);
  };

  bool runtimeDomainClampWarned = false;
  auto clampBodyBeforeMarkers = [&]() -> bool {
    const T halfLength = 0.5 * bodyLength;
    const T halfWidthWithWave = p.bodyRadius + tailAmpLU + T(2);
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

    if (xCm != oldX && (surgeOnlyMode || full3DofMode)) Vx = T(0);
    if (yCm != oldY && full3DofMode) Vy = T(0);
    if (xCm != oldX || yCm != oldY) {
      clampAudit.runtimeDomainClampHit = true;
      ++clampAudit.runtimeDomainClampCount;
      if (!runtimeDomainClampWarned) {
        clout << "WARNING: clamped body CM before marker generation from ("
              << oldX << ", " << oldY << ") to ("
              << xCm << ", " << yCm << ")." << std::endl;
        runtimeDomainClampWarned = true;
      }
      if (fixedPoseMode) {
        xCmRef = xCm;
        yCmRef = yCm;
      } else if (surgeOnlyMode && yCm != oldY) {
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
    const bool ok = (minX >= T(2) && maxX <= T(nx - 3) &&
                     minY >= T(2) && maxY <= T(ny - 3));
    if (!ok) {
      clout << "ERROR: IBM markers left the lattice during " << context
            << ". bounds x=[" << minX << ", " << maxX << "] y=["
            << minY << ", " << maxY << "] domain=[0," << (nx - 1)
            << "]x[0," << (ny - 1) << "]. Terminating cleanly."
            << std::endl;
    }
    return ok;
  };

  auto buildMarkersChecked = [&](T t, LagrangianMarkers& target,
                                 const char* context) -> bool {
    applyModeDefinitionState();
    if (!clampBodyBeforeMarkers()) return false;
    applyModeDefinitionState();
    buildLagrangianMarkers(
      p, t, Vx, Vy, omegaZ, xCm, yCm, theta, dtLbm,
      deformationAmpScale(), target);
    return markersInDomain(target, context);
  };

  // ---- Build initial markers ----
  LagrangianMarkers markers;
  if (!buildMarkersChecked(0.0, markers, "initialization")) {
    return 1;
  }

  // Kinematic time offset:
  //  * undulation: legacy behavior — gait clock during the main loop is
  //    offset by warmupDuration so that the gait is already past restTime
  //    and partially ramped when history starts.
  //  * rest / none: gait clock starts at 0 at history t=0, so history t=0
  //    coincides with the first real gait ramp onset.  This makes per-cycle
  //    averages from cycle 1 phase-aligned with the gait kinematics.
  const T warmupDuration =
    (warmupMode == WarmupMode::Undulation) ? T(p.nWarmup) * dtLbm : T(0);
  auto kinematicTimeAt = [&](T reportedTime) -> T {
    return warmupDuration + reportedTime;
  };

  // ---- Warmup ----
  IbmResult ibmRes;
  if (warmupMode == WarmupMode::None) {
    clout << "Warmup mode: none — skipping warmup loop." << std::endl;
  } else {
    const char* modeLabel =
      (warmupMode == WarmupMode::Rest) ? "rest" : "undulation";
    clout << "Starting warmup (" << p.nWarmup << " steps, mode="
          << modeLabel << ")..." << std::endl;
    for (int step = 0; step < p.nWarmup; ++step) {
      // For rest mode we always pass tWarmup=0: actuationProfile returns
      // ampScale=0 there, so the body stays in its still pose and the
      // fluid relaxes to a quiescent / boundary-consistent state.  For
      // undulation mode we advance the gait clock through the warmup so
      // the gait is partially ramped by the time history begins.
      const T tWarmup = (warmupMode == WarmupMode::Rest)
                      ? T(0) : T(step) * dtLbm;
      if (!buildMarkersChecked(tWarmup, markers, "warmup")) {
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
  std::vector<T> histT, histVx, histVy, histWz;
  std::vector<T> histFx, histFy, histTz;
  std::vector<T> histXcm, histYcm, histTheta;
  std::vector<T> histPower, histEta;
  // v6 power decomposition: P_total = P_rigid + P_def (per-frame, ds-weighted).
  std::vector<T> histPowerRigid, histPowerDef;
  // Marker-local residual slip after the multi-direct-forcing iterations.
  std::vector<T> histMeanResidualSlip, histMaxResidualSlip;
  // Body-frame diagnostics
  std::vector<T> histUswim, histUlat;
  // forwardNetForce: net hydrodynamic force projected onto the body forward
  // axis (NOT a separated propulsive thrust).  Renamed from the previous
  // misleading `Fthrust` / `histFthrust`.
  std::vector<T> histForwardNetForce, histFlat;
  std::vector<T> histRe, histSt, histMach, histUstar;
  // IBM diagnostics are kept separate from the rigid-body forces:
  //  mean/maxSlip          -> per-marker velocity-mismatch magnitudes
  //  mean/maxMarkerForce   -> per-marker line-force density magnitudes
  std::vector<T> histMeanSlip, histMaxSlip;
  std::vector<T> histMeanMarkerForce, histMaxMarkerForce, histNormalizedSlip;

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

  // Emergency Mach-safety caps for rigid-body motion seen by the IBM markers.
  // These are numerical protections, not part of the model definition:
  //  - translational cap:
  //      full3dof only -> limits sqrt(Vx^2 + Vy^2)
  //  - yaw-rate cap:
  //      full3dof only; surge_only defines omegaZ = 0 exactly
  //
  // The LBM Mach number should stay below ~0.1 (u < 0.058 lu/step).
  // 0.05 leaves a small safety margin without clipping ordinary low-speed runs.
  T maxBodyTranslationSpeedLat = 0.05;
  T maxBodyYawRateLat = 0.005;

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
    int ibmDiagSamples = 0;

    for (int sub = 0; sub < p.substeps; ++sub) {
      T tSub = tStart + sub * dtLbm;
      T tKinematics = kinematicTimeAt(tSub);

      // Coupling order: build the current geometry, evaluate IBM once,
      // collide-stream the fluid, then update the selected rigid-body DOFs.
      // Multiple IBM evaluations before one fluid update are avoided because
      // the interpolated fluid velocity would not change between evaluations.
      if (!buildMarkersChecked(tKinematics, markers, "main loop")) {
        return 1;
      }

      ibmStep(p, markers, xCm, yCm, Vx, Vy, omegaZ, sLattice,
              alphaIBM, ibmIterations, resetEulerForce, ibmRes);
      if (std::isfinite(ibmRes.meanSlipMag) &&
          std::isfinite(ibmRes.maxSlipMag) &&
          std::isfinite(ibmRes.meanMarkerForceMag) &&
          std::isfinite(ibmRes.maxMarkerForceMag) &&
          std::isfinite(ibmRes.meanDesiredMarkerSpeedMag)) {
        meanSlipAccum += ibmRes.meanSlipMag;
        meanMarkerForceAccum += ibmRes.meanMarkerForceMag;
        meanDesiredSpeedAccum += ibmRes.meanDesiredMarkerSpeedMag;
        maxSlipFrame = std::max(maxSlipFrame, ibmRes.maxSlipMag);
        maxMarkerForceFrame = std::max(maxMarkerForceFrame, ibmRes.maxMarkerForceMag);
        if (std::isfinite(ibmRes.meanResidualSlipMag) &&
            std::isfinite(ibmRes.maxResidualSlipMag)) {
          meanResidualSlipAccum += ibmRes.meanResidualSlipMag;
          maxResidualSlipFrame   = std::max(maxResidualSlipFrame,
                                            ibmRes.maxResidualSlipMag);
        }
        ++ibmDiagSamples;
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
        applyRigidBodyForceUpdate(simCase, tKinematics, p.restTime,
                                  fxS, fyS, tzS, bodyInertia, bodyState);
      }

      // Position/orientation: Δx = V, Δθ = ω  (δt = 1).
      advanceRigidBodyPose(simCase, {xCmRef, yCmRef, thetaRef}, bodyState);

      // In surge_only mode, we solve a reduced 1-DOF system:
      // Only Vx is integrated dynamically.
      // No artificial velocity clamping is applied,
      // to preserve physical consistency of force -> acceleration -> velocity.
      if (full3DofMode) {
        const T translationSpeedLat = std::sqrt(Vx * Vx + Vy * Vy);
        if (translationSpeedLat > maxBodyTranslationSpeedLat) {
          if (!clampAudit.speedClampHit) {
            clout << "WARNING: body-speed Mach safety cap activated at t=" << tSub
                  << "  |V|=" << translationSpeedLat
                  << " > " << maxBodyTranslationSpeedLat << std::endl;
            clampAudit.speedClampHit = true;
          }
          ++clampAudit.speedClampCount;
          const T factor = maxBodyTranslationSpeedLat / translationSpeedLat;
          Vx *= factor;
          Vy *= factor;
        }
        const T yawRateLat = std::abs(omegaZ);
        if (yawRateLat > maxBodyYawRateLat) {
          if (!clampAudit.omegaClampHit) {
            clout << "WARNING: yaw-rate safety cap activated at t=" << tSub
                  << "  |omega|=" << yawRateLat << " > " << maxBodyYawRateLat
                  << std::endl;
            clampAudit.omegaClampHit = true;
          }
          ++clampAudit.omegaClampCount;
          omegaZ *= maxBodyYawRateLat / yawRateLat;
        }
      }

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
    T diagnosticSpeed = fixedInflow ? std::abs(p.inflowVelocity) : absUswim;
    // Mach number (lattice).  fixed_inflow reports the imposed inflow speed.
    T machLat = (fixedInflow ? std::abs(p.inflowVelocity) : bodySpeed)
              / std::sqrt(T(1) / T(3));
    // Reynolds number: Re = U_ref * total capsule length / nu.
    T Re_inst = diagnosticSpeed * bodyLength / p.nu();
    // Strouhal number:  St = f * A_pp / |U_swim|
    //   f = eelFreq [Hz],  dtLbm [s/step],  A_pp [lu],  U [lu/step]
    //   St = (f * dtLbm) * A_pp / |U|   (all in lattice units)
    T St_inst = (!fixedInflow && absUswim > 1e-12)
              ? freqLat * tailPP / absUswim
              : 0.0;
    // Non-dimensional swimming speed:  U* = U_swim / (f * L)
    T Ustar = (!fixedInflow && absUswim > 1e-12 && freqLat > 1e-12)
            ? absUswim / (freqLat * bodyLength)
            : 0.0;

    // Store history
    const T tFrame = tStart + p.dtAnim;

    histT.push_back(tFrame);
    histVx.push_back(Vx);
    histVy.push_back(Vy);
    histWz.push_back(omegaZ);
    histFx.push_back(FxAvg);
    histFy.push_back(FyAvg);
    histTz.push_back(TzAvg);
    histXcm.push_back(xCm);
    histYcm.push_back(yCm);
    histTheta.push_back(theta);
    histPower.push_back(PAvg);
    histEta.push_back(eta);
    histUswim.push_back(Uswim);
    histUlat.push_back(Ulateral);
    histForwardNetForce.push_back(forwardNetForce);
    histFlat.push_back(Flat);
    histRe.push_back(Re_inst);
    histSt.push_back(St_inst);
    histMach.push_back(machLat);
    histUstar.push_back(Ustar);
    histMeanSlip.push_back(meanSlipFrame);
    histMaxSlip.push_back(maxSlipFrame);
    histMeanMarkerForce.push_back(meanMarkerForceFrame);
    histMaxMarkerForce.push_back(maxMarkerForceFrame);
    histNormalizedSlip.push_back(normalizedSlipFrame);
    histPowerRigid.push_back(PRigidAvg);
    histPowerDef.push_back(PDefAvg);
    histMeanResidualSlip.push_back(meanResidualSlipFrame);
    histMaxResidualSlip.push_back(maxResidualSlipFrame);

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

    // Ensure exported fields are synchronized with the latest streamed state.
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Export gating: respect exportInterval and run mode
    const bool doExport = writeSpatialExports && (frame % exportInterval == 0);
    const bool doDiagnostics = doExport && (runMode == RunMode::Full);

    if (doExport) {
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

      // Compute body mask once, share across all VTI writers
      auto bodyMask = computeBodyMask(frameMarkers, nx, ny);
      exportFields.sample(sLattice, bodyMask, doDiagnostics);

      writeVelocityVTI(velDir, frame, vtiWriter, exportFields);
      writeVorticityVTI(vortDir, frame, vtiWriter, exportFields);
      if (doDiagnostics) {
        writeDiagnosticsVTI(diagDir, frame, vtiWriter, exportFields);
      }

      if (doDiagnostics) resetEulerForce();

      writeBodySnapshotCSV(bodyOutDir, frame, frameMarkers);
      writeBodyVTP(bodyOutDir, frame, p, frameMarkers);
    }
  }

  // ============================================================
  //  Save ParaView .pvd time-series files (one per field group)
  // ============================================================
  if (writeSpatialExports) {
    std::string pvdPath = writePVD(runOutDir, "eel3dof_velocity.pvd", "vtkData/velocity/", "eel3dof_velocity_", ".vti", nFrames, exportInterval, p.dtAnim);
    clout << "Saved PVD: " << pvdPath << std::endl;
    pvdPath = writePVD(runOutDir, "eel3dof_vorticity.pvd", "vtkData/vorticity/", "eel3dof_vorticity_", ".vti", nFrames, exportInterval, p.dtAnim);
    clout << "Saved PVD: " << pvdPath << std::endl;
    if (runMode == RunMode::Full) {
      pvdPath = writePVD(runOutDir, "eel3dof_diagnostics.pvd", "vtkData/diagnostics/", "eel3dof_diagnostics_", ".vti", nFrames, exportInterval, p.dtAnim);
      clout << "Saved PVD: " << pvdPath << std::endl;
    }
    pvdPath = writePVD(runOutDir, "eel3dof_body.pvd", "bodyData/", "eel3dof_body_", ".vtp", nFrames, exportInterval, p.dtAnim);
    clout << "Saved PVD: " << pvdPath << std::endl;
  } else {
    clout << "Verification mode: skipped VTK/body snapshot exports." << std::endl;
  }

  // ============================================================
  //  Append one-run aspect-ratio summary
  // ============================================================
  SteadySummary steady = computeSteadySummary(
    histT, histUswim, histPower, histForwardNetForce, histFlat,
    histRe, histSt, histUstar, histEta,
    histMeanSlip, histMaxSlip, histMeanMarkerForce, histMaxMarkerForce,
    histNormalizedSlip,
    histPowerRigid, histPowerDef,
    histMeanResidualSlip, histMaxResidualSlip,
    tCut, mass);

  if (steady.nSamples == 0) {
    clout << "WARNING: no steady-state samples found after tCut=" << tCut
          << ". Increase Ttotal or lower --tCut for summary metrics." << std::endl;
  }

  const int recommendedSteadySamples = recommendedSteadySampleCount(p, simCase);
  const bool enoughSteadySamples = (steady.nSamples >= recommendedSteadySamples);

  // ----------------------------------------------------------------
  //  Cycle-averaged steady-state diagnostics
  // ----------------------------------------------------------------
  //  period = 1 / eelFreq (physical seconds, same units as histT).  Bins
  //  start at the first sample with t > tCut.  Convergence is reported over
  //  the last min(5, nCycles) complete cycles.
  const T cyclePeriod = (p.eelFreq > T(1e-12)) ? (T(1) / p.eelFreq) : T(0);
  std::vector<CycleAverage> cycles = computeCycleAverages(
    histT, histUswim, histPower, histForwardNetForce,
    histRe, histSt, histUstar,
    histMeanSlip, histNormalizedSlip,
    tCut, cyclePeriod, mass);
  const CycleConvergence cycleConv = computeCycleConvergence(cycles);

  SummaryCsvInputs summaryInput;
  summaryInput.p = p;
  summaryInput.steady = steady;
  summaryInput.cycleConv = cycleConv;
  summaryInput.simCase = simCase;
  summaryInput.warmupMode = warmupMode;
  summaryInput.gaitNormalization = effectiveGaitNormalization;
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
  summaryInput.effectiveEelFreq = effectiveEelFreq;
  summaryInput.effectiveEelA0 = effectiveEelA0;
  summaryInput.alphaIBM = alphaIBM;
  summaryInput.legacyKappaInputUsed = legacyKappaInputUsed;
  summaryInput.ibmIterations = ibmIterations;

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
  verificationInput.gaitNormalization = effectiveGaitNormalization;
  verificationInput.initialPlacementClamped = clampAudit.initialPlacementClamped;
  verificationInput.initialPlacementClampCount = clampAudit.initialPlacementClampCount;
  verificationInput.speedClampHit = clampAudit.speedClampHit;
  verificationInput.runtimeDomainClampHit = clampAudit.runtimeDomainClampHit;
  verificationInput.speedClampCount = clampAudit.speedClampCount;
  verificationInput.runtimeDomainClampCount = clampAudit.runtimeDomainClampCount;
  verificationInput.omegaClampCount = clampAudit.omegaClampCount;
  verificationInput.nx = nx;
  verificationInput.ny = ny;
  verificationInput.dtLbm = dtLbm;
  verificationInput.tCut = tCut;
  verificationInput.bodyCenterlineLength = bodyCenterlineLength;
  verificationInput.bodyLength = bodyLength;
  verificationInput.bodyWidth = bodyWidth;
  verificationInput.mass = mass;
  verificationInput.Ibody = Ibody;
  verificationInput.effectiveEelFreq = effectiveEelFreq;
  verificationInput.effectiveEelA0 = effectiveEelA0;
  verificationInput.alphaIBM = alphaIBM;
  verificationInput.legacyKappaInputUsed = legacyKappaInputUsed;
  verificationInput.ibmIterations = ibmIterations;

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
        << "  kappa=" << p.kappa
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
  clout << "speedClampCount=" << clampAudit.speedClampCount
        << "  runtimeDomainClampHit=" << (clampAudit.runtimeDomainClampHit ? 1 : 0)
        << "  runtimeDomainClampCount=" << clampAudit.runtimeDomainClampCount
        << "  omegaClampCount=" << clampAudit.omegaClampCount << std::endl;
  clout << "initialPlacementClamped=" << (clampAudit.initialPlacementClamped ? 1 : 0)
        << "  initialPlacementClampCount=" << clampAudit.initialPlacementClampCount
        << "  (one-shot setup nudge from initialPositionFactor / rightMargin; "
        << "not a runtime issue)" << std::endl;
  clout << "warmupMode=" << warmupModeName(warmupMode)
        << "  gaitNormalization=" << gaitNormalizationName(effectiveGaitNormalization)
        << "  effectiveEelFreq=" << effectiveEelFreq
        << "  effectiveEelA0=" << effectiveEelA0 << std::endl;
  clout << "IBM(v6): alphaIBM=" << alphaIBM
        << "  ibmIterations=" << ibmIterations
        << "  legacyKappaInput=" << (legacyKappaInputUsed ? 1 : 0)
        << "  meanResidualSlip=" << steady.meanResidualSlip
        << "  maxResidualSlip=" << steady.maxResidualSlip << std::endl;
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
    HistoryCsvData historyCsv;
    historyCsv.histT = &histT;
    historyCsv.histVx = &histVx;
    historyCsv.histVy = &histVy;
    historyCsv.histWz = &histWz;
    historyCsv.histFx = &histFx;
    historyCsv.histFy = &histFy;
    historyCsv.histTz = &histTz;
    historyCsv.histXcm = &histXcm;
    historyCsv.histYcm = &histYcm;
    historyCsv.histTheta = &histTheta;
    historyCsv.histPower = &histPower;
    historyCsv.histEta = &histEta;
    historyCsv.histUswim = &histUswim;
    historyCsv.histUlat = &histUlat;
    historyCsv.histForwardNetForce = &histForwardNetForce;
    historyCsv.histFlat = &histFlat;
    historyCsv.histRe = &histRe;
    historyCsv.histSt = &histSt;
    historyCsv.histMach = &histMach;
    historyCsv.histUstar = &histUstar;
    historyCsv.histMeanSlip = &histMeanSlip;
    historyCsv.histMaxSlip = &histMaxSlip;
    historyCsv.histMeanMarkerForce = &histMeanMarkerForce;
    historyCsv.histMaxMarkerForce = &histMaxMarkerForce;
    historyCsv.histNormalizedSlip = &histNormalizedSlip;
    historyCsv.histPowerRigid = &histPowerRigid;
    historyCsv.histPowerDef = &histPowerDef;
    historyCsv.histMeanResidualSlip = &histMeanResidualSlip;
    historyCsv.histMaxResidualSlip = &histMaxResidualSlip;
    const std::string fname = writeHistoryCsv(runOutDir, historyCsv);
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
