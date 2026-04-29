#pragma once

#include <string>

// Run modes - control which outputs are written and at what cadence.
enum class RunMode { Preview, Standard, Full };

inline RunMode parseRunMode(const std::string& s) {
  if (s == "preview")  return RunMode::Preview;
  if (s == "standard") return RunMode::Standard;
  return RunMode::Full;
}

// Study modes - control verification bookkeeping/output around unchanged physics.
enum class StudyMode { Standard, Verification };

inline StudyMode parseStudyMode(const std::string& s) {
  if (s == "verification") return StudyMode::Verification;
  return StudyMode::Standard;
}

inline const char* studyModeName(StudyMode mode) {
  return (mode == StudyMode::Verification) ? "verification" : "standard";
}

// Simulation cases - control physical degrees of freedom.
enum class SimulationCase { FixedInflow, FixedUndulation, SurgeOnly, Full3Dof };

inline SimulationCase parseSimulationCase(const std::string& s) {
  if (s == "fixed_inflow")      return SimulationCase::FixedInflow;
  if (s == "fixed_undulation")  return SimulationCase::FixedUndulation;
  if (s == "full3dof")          return SimulationCase::Full3Dof;
  return SimulationCase::SurgeOnly;
}

inline const char* simulationCaseName(SimulationCase simCase) {
  switch (simCase) {
    case SimulationCase::FixedInflow:     return "fixed_inflow";
    case SimulationCase::FixedUndulation: return "fixed_undulation";
    case SimulationCase::SurgeOnly:       return "surge_only";
    case SimulationCase::Full3Dof:        return "full3dof";
  }
  return "surge_only";
}

inline bool parseBool(const std::string& s) {
  return s == "1" || s == "true" || s == "yes" || s == "on";
}

// Warmup modes.
enum class WarmupMode { Rest, Undulation, None };

inline WarmupMode parseWarmupMode(const std::string& s) {
  if (s == "rest")       return WarmupMode::Rest;
  if (s == "undulation") return WarmupMode::Undulation;
  if (s == "none")       return WarmupMode::None;
  return WarmupMode::Rest;
}

inline const char* warmupModeName(WarmupMode m) {
  switch (m) {
    case WarmupMode::Rest:       return "rest";
    case WarmupMode::Undulation: return "undulation";
    case WarmupMode::None:       return "none";
  }
  return "rest";
}

// Gait normalization modes.
enum class GaitNormalization { Fixed, TailAmpRatio, TargetSt };

inline GaitNormalization parseGaitNormalization(const std::string& s) {
  if (s == "tailAmpRatio") return GaitNormalization::TailAmpRatio;
  if (s == "targetSt")     return GaitNormalization::TargetSt;
  return GaitNormalization::Fixed;
}

inline const char* gaitNormalizationName(GaitNormalization g) {
  switch (g) {
    case GaitNormalization::Fixed:        return "fixed";
    case GaitNormalization::TailAmpRatio: return "tailAmpRatio";
    case GaitNormalization::TargetSt:     return "targetSt";
  }
  return "fixed";
}
