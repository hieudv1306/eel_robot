#include "physics/diagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <utility>



// ============================================================
//  Cycle-averaged steady-state diagnostics
// ============================================================
//  A "cycle" is one undulation period, period = 1 / eelFreq, expressed in the
//  same physical-time coordinate as histT.  Cycle bins start at the first
//  sample with t > tCut and march forward in equal-period intervals.  Only
//  bins that lie entirely inside the recorded history are kept.  Per-cycle
//  averages are arithmetic means of finite samples falling in each bin;
//  derived ratios (CoT, hydroCost, transportEfficiencyProxy) are computed
//  from the per-cycle averaged Uswim and power, not from per-frame ratios.
//
//  Assumptions used in cycle averaging:
//    * The first cycle bin starts at the first frame with t > tCut, so its
//      starting phase is whatever phase the body happens to have at tCut.
//      For long enough runs this does not bias the mean.
//    * Only complete bins are kept (a bin is dropped if its [tStart, tEnd]
//      window extends past the last recorded frame).
//    * Cycle-averaged power uses |power| (input cost), matching the steady
//      window definition.




// Minimal adequacy check for steady-window averaging in verification sweeps.
// This is only an audit flag; it is not a substitute for convergence testing.
int recommendedSteadySampleCount(const EelParams& p, SimulationCase simCase)
{
  if (simCase == SimulationCase::FixedInflow) {
    return 10;
  }
  if (!(p.dtAnim > T(0))) {
    return 10;
  }
  if (p.eelFreq > T(1e-12)) {
    const T framesPerPeriod = T(1) / (p.eelFreq * p.dtAnim);
    return std::max(10, static_cast<int>(std::ceil(T(2) * framesPerPeriod)));
  }
  return 10;
}

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
    T mass)
{
  SteadySummary out;
  const T eps = T(1e-12);
  // v6 accumulators.  These are kept as locals so SteadySummary stays a POD
  // of finalized averages (matching v5 conventions).
  T sumConstraintPower = T(0);
  T sumRigidBodyPower  = T(0);
  T sumDeformationPower = T(0);
  T sumAbsConstraintPower = T(0);
  T sumAbsRigidBodyPower  = T(0);
  T sumAbsDeformationPower = T(0);
  T sumMeanResidual = T(0);
  T maxResidualWindow = T(0);
  int residualSamples = 0;
  int powerSamples = 0;

  for (size_t i = 0; i < t.size(); ++i) {
    if (!(t[i] > tCut)) continue;
    if (i >= meanSlip.size() || i >= maxSlip.size() ||
        i >= meanMarkerForce.size() || i >= maxMarkerForce.size() ||
        i >= normalizedSlip.size()) {
      continue;
    }
    const T u  = uswim[i];
    const T p  = std::abs(power[i]);  // input power is reported as positive cost
    // Net hydrodynamic force on the body projected onto the forward swim
    // axis.  In surge_only at theta=0 this is -Fx.  This is NOT a separated
    // propulsive thrust component.
    const T fNet = forwardNetForce[i];
    const T fl = lateralForce[i];
    if (!std::isfinite(u) || !std::isfinite(p) || !std::isfinite(fNet) ||
        !std::isfinite(fl) || !std::isfinite(re[i]) || !std::isfinite(st[i]) ||
        !std::isfinite(ustar[i]) || !std::isfinite(eta[i]) ||
        !std::isfinite(meanSlip[i]) || !std::isfinite(maxSlip[i]) ||
        !std::isfinite(meanMarkerForce[i]) || !std::isfinite(maxMarkerForce[i]) ||
        !std::isfinite(normalizedSlip[i])) {
      continue;
    }

    ++out.nSamples;
    out.meanU += u;
    out.meanP += p;
    out.meanForwardNetForce += fNet;
    out.meanAbsForwardNetForce += std::abs(fNet);
    out.meanSignedForwardNetForce += fNet;
    out.meanLateralForce += fl;
    out.meanRe += re[i];
    out.meanSt += st[i];
    out.meanUstar += ustar[i];
    out.meanEta += eta[i];
    out.meanSlip += meanSlip[i];
    out.meanMaxSlip += maxSlip[i];
    out.meanMarkerForce += meanMarkerForce[i];
    out.meanMaxMarkerForce += maxMarkerForce[i];
    out.meanNormalizedSlip += normalizedSlip[i];

    // v6: power decomposition averaged from per-frame ds-weighted sums.
    if (i < powerRigid.size() && i < powerDef.size()) {
      const T pTot = power[i];           // signed (NOT abs)
      const T pRig = powerRigid[i];
      const T pDef = powerDef[i];
      if (std::isfinite(pTot) && std::isfinite(pRig) && std::isfinite(pDef)) {
        sumConstraintPower    += pTot;
        sumRigidBodyPower     += pRig;
        sumDeformationPower   += pDef;
        sumAbsConstraintPower += std::abs(pTot);
        sumAbsRigidBodyPower  += std::abs(pRig);
        sumAbsDeformationPower += std::abs(pDef);
        ++powerSamples;
      }
    }
    if (i < meanResidualSlipPerFrame.size() &&
        i < maxResidualSlipPerFrame.size()) {
      const T mr = meanResidualSlipPerFrame[i];
      const T xr = maxResidualSlipPerFrame[i];
      if (std::isfinite(mr) && std::isfinite(xr)) {
        sumMeanResidual += mr;
        maxResidualWindow = std::max(maxResidualWindow, xr);
        ++residualSamples;
      }
    }
  }

  if (out.nSamples == 0) {
    return out;
  }

  const T invN = T(1) / T(out.nSamples);
  out.meanU *= invN;
  out.meanP *= invN;
  out.meanForwardNetForce *= invN;
  out.meanLateralForce *= invN;
  out.meanRe *= invN;
  out.meanSt *= invN;
  out.meanUstar *= invN;
  out.meanEta *= invN;
  out.meanSlip *= invN;
  out.meanMaxSlip *= invN;
  out.meanMarkerForce *= invN;
  out.meanMaxMarkerForce *= invN;
  out.meanNormalizedSlip *= invN;
  out.meanAbsForwardNetForce *= invN;
  out.meanSignedForwardNetForce *= invN;

  // meanSignedForwardNetForce is the same quantity as meanForwardNetForce
  // (signed net force projected onto the forward axis).  The pair is exposed
  // explicitly in the output schema for clarity against meanAbsForwardNetForce.
  out.meanSignedForwardNetForce = out.meanForwardNetForce;

  // etaNetForceDiagnostic is only a diagnostic of work done against the mean
  // net forward force.  For self-propelled steady swimming the mean net
  // forward force may approach zero or change sign, so this ratio is NOT a
  // reliable propulsion efficiency.  Use CoT, hydroCost,
  // transportEfficiencyProxy, and meanU for aspect-ratio ranking instead.
  if (out.meanP > eps) {
    out.etaNetForceDiagnostic = (out.meanForwardNetForce * out.meanU) / out.meanP;
    out.transportEfficiencyProxy = out.meanU / out.meanP;
  }
  if (mass > eps && out.meanU > eps) {
    out.CoT = out.meanP / (mass * out.meanU);
  }
  if (out.meanU > eps) {
    out.hydroCost = out.meanP / out.meanU;
  }

  // v6 — finalize the power decomposition over the steady window.  All three
  // averages share the same sample base so the additive identity
  // P_total ≈ P_rigid + P_def is preserved to floating-point round-off.
  if (powerSamples > 0) {
    const T invP = T(1) / T(powerSamples);
    out.meanConstraintPower    = sumConstraintPower    * invP;
    out.meanRigidBodyPower     = sumRigidBodyPower     * invP;
    out.meanDeformationPower   = sumDeformationPower   * invP;
    out.meanAbsConstraintPower = sumAbsConstraintPower * invP;
    out.meanAbsRigidBodyPower  = sumAbsRigidBodyPower  * invP;
    out.meanAbsDeformationPower = sumAbsDeformationPower * invP;
  }
  if (residualSamples > 0) {
    out.meanResidualSlip = sumMeanResidual / T(residualSamples);
  }
  out.maxResidualSlip = maxResidualWindow;

  // Deformation-only efficiency proxies (preferred metric for AR ranking
  // because they isolate gait energetics from translational kinetic-energy
  // bookkeeping).
  if (out.meanAbsDeformationPower > eps) {
    out.transportEfficiencyDef = out.meanU / out.meanAbsDeformationPower;
    out.cotDef = out.meanAbsDeformationPower / std::max(out.meanU, eps);
  }

  return out;
}

// Build per-cycle averages over the steady window [tCut, t.back()] using
// undulation period = 1 / eelFreq.  See CycleAverage docstring for the binning
// and assumption notes.  Returns an empty vector if period <= 0 or there is
// not even one complete cycle inside the steady window.
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
    T mass)
{
  std::vector<CycleAverage> cycles;
  if (!(period > T(0)) || t.empty()) {
    return cycles;
  }

  size_t i0 = 0;
  while (i0 < t.size() && !(t[i0] > tCut)) ++i0;
  if (i0 >= t.size()) return cycles;

  const T tStart0 = t[i0];
  const T tLast   = t.back();
  const int nCycles = static_cast<int>(std::floor((tLast - tStart0) / period));
  if (nCycles <= 0) return cycles;

  cycles.reserve(nCycles);
  const T eps = T(1e-12);

  size_t cursor = i0;
  for (int c = 0; c < nCycles; ++c) {
    const T cStart = tStart0 + T(c) * period;
    const T cEnd   = tStart0 + T(c + 1) * period;

    CycleAverage ca;
    ca.tStart = cStart;
    ca.tEnd   = cEnd;

    int n = 0;
    T sumU = 0, sumP = 0, sumFnet = 0, sumRe = 0, sumSt = 0, sumUstar = 0;
    T sumSlip = 0, sumNormSlip = 0;

    size_t i = cursor;
    for (; i < t.size(); ++i) {
      const T ti = t[i];
      if (ti <= cStart) continue;
      if (ti > cEnd) break;
      if (i >= meanSlip.size() || i >= normalizedSlip.size()) continue;
      const T u = uswim[i];
      const T pwr = std::abs(power[i]);
      const T fNet = forwardNetForce[i];
      if (!std::isfinite(u) || !std::isfinite(pwr) || !std::isfinite(fNet) ||
          !std::isfinite(re[i]) || !std::isfinite(st[i]) ||
          !std::isfinite(ustar[i]) || !std::isfinite(meanSlip[i]) ||
          !std::isfinite(normalizedSlip[i])) {
        continue;
      }
      ++n;
      sumU += u;
      sumP += pwr;
      sumFnet += fNet;
      sumRe += re[i];
      sumSt += st[i];
      sumUstar += ustar[i];
      sumSlip += meanSlip[i];
      sumNormSlip += normalizedSlip[i];
    }
    cursor = i;

    if (n == 0) {
      continue;
    }

    const T inv = T(1) / T(n);
    ca.nSamples = n;
    ca.Uswim = sumU * inv;
    ca.power = sumP * inv;
    ca.forwardNetForce = sumFnet * inv;
    ca.Re = sumRe * inv;
    ca.St = sumSt * inv;
    ca.Ustar = sumUstar * inv;
    ca.meanSlip = sumSlip * inv;
    ca.normalizedSlip = sumNormSlip * inv;

    if (mass > eps && ca.Uswim > eps) {
      ca.CoT = ca.power / (mass * ca.Uswim);
    }
    if (ca.Uswim > eps) {
      ca.hydroCost = ca.power / ca.Uswim;
    }
    if (ca.power > eps) {
      ca.transportEfficiencyProxy = ca.Uswim / ca.power;
    }

    cycles.push_back(ca);
  }

  return cycles;
}

// Compute mean and coefficient of variation of selected per-cycle quantities
// over the last min(5, cycles.size()) cycles.  cycleConverged is true only if
// at least 3 cycles are available AND the CoV of meanUstar, CoT, and
// hydroCost are each below 0.05.
CycleConvergence computeCycleConvergence(const std::vector<CycleAverage>& cycles)
{
  CycleConvergence out;
  const int total = static_cast<int>(cycles.size());
  if (total <= 0) return out;

  const int nUse = std::min(5, total);
  out.nSteadyCycles = nUse;
  const int start = total - nUse;

  auto meanCv = [&](auto extractor) -> std::pair<T, T> {
    T sum = 0;
    int count = 0;
    for (int i = start; i < total; ++i) {
      const T v = extractor(cycles[i]);
      if (!std::isfinite(v)) continue;
      sum += v;
      ++count;
    }
    if (count == 0) return {T(0), T(0)};
    const T mean = sum / T(count);
    T var = 0;
    for (int i = start; i < total; ++i) {
      const T v = extractor(cycles[i]);
      if (!std::isfinite(v)) continue;
      const T d = v - mean;
      var += d * d;
    }
    var /= T(count);
    const T sd = std::sqrt(var);
    const T cv = (std::abs(mean) > T(1e-12)) ? sd / std::abs(mean) : T(0);
    return {mean, cv};
  };

  auto [muUstar, cvUstar] = meanCv([](const CycleAverage& c) { return c.Ustar; });
  auto [muCoT,   cvCoT]   = meanCv([](const CycleAverage& c) { return c.CoT; });
  auto [muHydro, cvHydro] = meanCv([](const CycleAverage& c) { return c.hydroCost; });
  auto [muNslip, cvNslip] = meanCv([](const CycleAverage& c) { return c.normalizedSlip; });

  out.cycleMeanUstar = muUstar;
  out.cycleCvUstar = cvUstar;
  out.cycleMeanCoT = muCoT;
  out.cycleCvCoT = cvCoT;
  out.cycleMeanHydroCost = muHydro;
  out.cycleCvHydroCost = cvHydro;
  out.cycleMeanNormalizedSlip = muNslip;
  out.cycleCvNormalizedSlip = cvNslip;

  out.cycleConverged = (out.nSteadyCycles >= 3 &&
                        out.cycleCvUstar < T(0.05) &&
                        out.cycleCvCoT < T(0.05) &&
                        out.cycleCvHydroCost < T(0.05));
  return out;
}

