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

// Direction of the travelling actuation wave along the spine parameter s.
// The spine uses s=0 at the head and s=L at the tail.
enum class WaveDirection { HeadToTail, TailToHead };

inline WaveDirection parseWaveDirection(const std::string& s) {
  if (s == "tail_to_head" || s == "tailToHead" ||
      s == "tail-head" || s == "reverse" || s == "-1") {
    return WaveDirection::TailToHead;
  }
  return WaveDirection::HeadToTail;
}

inline const char* waveDirectionName(WaveDirection d) {
  return (d == WaveDirection::TailToHead)
       ? "tail_to_head" : "head_to_tail";
}

inline double waveDirectionSign(WaveDirection d) {
  return (d == WaveDirection::TailToHead) ? -1.0 : 1.0;
}

// Marker/body deformation source.
// PrescribedWave keeps the legacy direct gait-to-marker path.
// SoftBackbone builds markers from a soft-backbone state that follows the
// preferred-curvature wave; fluid-to-backbone dynamics are added in the next
// coupling step.
enum class BodyKinematics { PrescribedWave, SoftBackbone };

inline BodyKinematics parseBodyKinematics(const std::string& s) {
  if (s == "soft_backbone" || s == "softBackbone" ||
      s == "backbone" || s == "soft") {
    return BodyKinematics::SoftBackbone;
  }
  return BodyKinematics::PrescribedWave;
}

inline const char* bodyKinematicsName(BodyKinematics k) {
  return (k == BodyKinematics::SoftBackbone)
       ? "soft_backbone" : "prescribed_wave";
}

// Soft-backbone time-integration scheme.
//   Overdamped : 1st-order relaxation toward the (fluid-offset) preferred
//                joint angles.  Cheap, robust, but ignores segment inertia
//                and only sees fluid loads as a quasi-static curvature shift.
//   Implicit   : 2nd-order Newton-Euler on segment angles with implicit
//                Euler in time; tridiagonal solve, unconditionally stable
//                for the stiff Dragon Skin K_theta.  Captures inertia and
//                phase lag with respect to the gait.
enum class SoftBackboneIntegrator { Overdamped, Implicit };

inline SoftBackboneIntegrator parseSoftBackboneIntegrator(const std::string& s) {
  if (s == "overdamped" || s == "first_order" || s == "1") {
    return SoftBackboneIntegrator::Overdamped;
  }
  return SoftBackboneIntegrator::Implicit;
}

inline const char* softBackboneIntegratorName(SoftBackboneIntegrator s) {
  return (s == SoftBackboneIntegrator::Overdamped)
       ? "overdamped" : "implicit";
}

// Body material used for inertia and soft-body stiffness estimates.
enum class BodyMaterial { NeutralLattice, DragonSkin20 };

inline BodyMaterial parseBodyMaterial(const std::string& s) {
  if (s == "neutral" || s == "neutral_lattice" || s == "lattice") {
    return BodyMaterial::NeutralLattice;
  }
  if (s == "dragon_skin_20" || s == "dragonskin20" ||
      s == "dragonSkin20" || s == "ds20") {
    return BodyMaterial::DragonSkin20;
  }
  return BodyMaterial::DragonSkin20;
}

inline const char* bodyMaterialName(BodyMaterial m) {
  return (m == BodyMaterial::DragonSkin20)
       ? "dragon_skin_20" : "neutral_lattice";
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
