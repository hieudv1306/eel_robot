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

// Top/bottom wall boundary type. NoSlip = bounce-back (legacy default);
// FreeSlip = full-slip (specular reflection) — removes channel-blockage drag
// from finite-ny domains, important for fair shape comparisons in AR sweeps.
enum class WallBoundary { NoSlip, FreeSlip };

inline WallBoundary parseWallBoundary(const std::string& s) {
  if (s == "freeslip" || s == "free_slip" || s == "fullslip") return WallBoundary::FreeSlip;
  return WallBoundary::NoSlip;
}

inline const char* wallBoundaryName(WallBoundary w) {
  return (w == WallBoundary::FreeSlip) ? "freeslip" : "noslip";
}

// Centerline kinematic construction.
// HeightWave preserves the legacy y=A(s)sin(...) geometry.
// InextensibleWave integrates a tangent-angle wave so segment lengths stay
// fixed while retaining the same wave phase/envelope controls.
enum class GeometryKinematics { HeightWave, InextensibleWave };

inline GeometryKinematics parseGeometryKinematics(const std::string& s) {
  if (s == "inextensible_wave" || s == "inextensible" ||
      s == "tangent_angle" || s == "tangentAngle") {
    return GeometryKinematics::InextensibleWave;
  }
  return GeometryKinematics::HeightWave;
}

inline const char* geometryKinematicsName(GeometryKinematics g) {
  return (g == GeometryKinematics::InextensibleWave)
       ? "inextensible_wave" : "height_wave";
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
