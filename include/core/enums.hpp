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
  if (s == "verification" || s == "ar_sweep" || s == "sweep") {
    return StudyMode::Verification;
  }
  return StudyMode::Standard;
}

inline const char* studyModeName(StudyMode mode) {
  return (mode == StudyMode::Verification) ? "verification" : "standard";
}

// Simulation cases.  Reduced to the surge-only production case; verification
// modes (fixed_inflow / fixed_undulation / full3dof) were retired with the
// soft-backbone-only research scope.
enum class SimulationCase { SurgeOnly };

inline SimulationCase parseSimulationCase(const std::string& /*s*/) {
  return SimulationCase::SurgeOnly;
}

inline const char* simulationCaseName(SimulationCase /*simCase*/) {
  return "surge_only";
}

inline bool parseBool(const std::string& s) {
  return s == "1" || s == "true" || s == "yes" || s == "on";
}

// Warmup modes.  "rest" lets the body sit still while the fluid relaxes;
// "none" skips the warmup loop entirely.
enum class WarmupMode { Rest, None };

inline WarmupMode parseWarmupMode(const std::string& s) {
  if (s == "none") return WarmupMode::None;
  return WarmupMode::Rest;
}

inline const char* warmupModeName(WarmupMode m) {
  return (m == WarmupMode::None) ? "none" : "rest";
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
