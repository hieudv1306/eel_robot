#pragma once

#include "core/enums.hpp"
#include "core/params.hpp"
#include "core/types.hpp"

#include <vector>

struct SteadySummary {
  int nSamples = 0;
  T meanU = 0.0;
  T meanP = 0.0;
  // Mean of the net hydrodynamic force projected onto the forward swim axis.
  // This is NOT a separated propulsive thrust component — it includes drag,
  // pressure, and viscous reactions and can vanish or change sign in steady
  // self-propelled swimming.  Use CoT / hydroCost / transportEfficiencyProxy
  // for propulsion ranking; keep meanForwardNetForce only as a force-balance
  // diagnostic.
  T meanForwardNetForce = 0.0;
  T meanLateralForce = 0.0;
  T meanRe = 0.0;
  T meanSt = 0.0;
  T meanUstar = 0.0;
  T meanEta = 0.0;
  T meanSlip = 0.0;
  T meanMaxSlip = 0.0;
  T meanMarkerForce = 0.0;
  T meanMaxMarkerForce = 0.0;
  T meanNormalizedSlip = 0.0;
  // Diagnostic-only: (meanForwardNetForce * meanU) / meanP.  In self-propelled
  // steady swimming the mean net force can approach zero, so this ratio is
  // unreliable as a propulsion efficiency.  Reported under the renamed field
  // etaNetForceDiagnostic for clarity.
  T etaNetForceDiagnostic = 0.0;
  T CoT = 0.0;
  T hydroCost = 0.0;
  T transportEfficiencyProxy = 0.0;
  T meanAbsForwardNetForce = 0.0;
  T meanSignedForwardNetForce = 0.0;
  // v6 power decomposition (signed and abs averages over the steady window).
  T meanConstraintPower = 0.0;
  T meanRigidBodyPower = 0.0;
  T meanDeformationPower = 0.0;
  T meanAbsConstraintPower = 0.0;
  T meanAbsRigidBodyPower = 0.0;
  T meanAbsDeformationPower = 0.0;
  // Marker-local residual slip from MDF iterations.
  T meanResidualSlip = 0.0;
  T maxResidualSlip  = 0.0;
  // Deformation-only efficiency proxies (preferred for AR ranking).
  T transportEfficiencyDef = 0.0;
  T cotDef = 0.0;
};

struct CycleAverage {
  T tStart = 0.0;
  T tEnd   = 0.0;
  int nSamples = 0;
  T Uswim = 0.0;
  T power = 0.0;
  T forwardNetForce = 0.0;
  T Re = 0.0;
  T St = 0.0;
  T Ustar = 0.0;
  T CoT = 0.0;
  T hydroCost = 0.0;
  T transportEfficiencyProxy = 0.0;
  T meanSlip = 0.0;
  T normalizedSlip = 0.0;
};

struct CycleConvergence {
  int nSteadyCycles = 0;
  T cycleMeanUstar = 0.0;
  T cycleCvUstar = 0.0;
  T cycleMeanCoT = 0.0;
  T cycleCvCoT = 0.0;
  T cycleMeanHydroCost = 0.0;
  T cycleCvHydroCost = 0.0;
  T cycleMeanNormalizedSlip = 0.0;
  T cycleCvNormalizedSlip = 0.0;
  bool cycleConverged = false;
};

int recommendedSteadySampleCount(const EelParams& p, SimulationCase simCase);

SteadySummary computeSteadySummary(
    const std::vector<T>& t,
    const std::vector<T>& uswim,
    const std::vector<T>& power,
    const std::vector<T>& forwardNetForce,
    const std::vector<T>& lateralForce,
    const std::vector<T>& re,
    const std::vector<T>& st,
    const std::vector<T>& ustar,
    const std::vector<T>& eta,
    const std::vector<T>& meanSlip,
    const std::vector<T>& maxSlip,
    const std::vector<T>& meanMarkerForce,
    const std::vector<T>& maxMarkerForce,
    const std::vector<T>& normalizedSlip,
    const std::vector<T>& powerRigid,
    const std::vector<T>& powerDef,
    const std::vector<T>& meanResidualSlipPerFrame,
    const std::vector<T>& maxResidualSlipPerFrame,
    T tCut,
    T mass);

std::vector<CycleAverage> computeCycleAverages(
    const std::vector<T>& t,
    const std::vector<T>& uswim,
    const std::vector<T>& power,
    const std::vector<T>& forwardNetForce,
    const std::vector<T>& re,
    const std::vector<T>& st,
    const std::vector<T>& ustar,
    const std::vector<T>& meanSlip,
    const std::vector<T>& normalizedSlip,
    T tCut,
    T period,
    T mass);

CycleConvergence computeCycleConvergence(const std::vector<CycleAverage>& cycles);
