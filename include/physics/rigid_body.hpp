#pragma once

#include "core/enums.hpp"
#include "core/params.hpp"
#include "core/types.hpp"

#include <cstdint>

struct BodyState {
  T xCm = 0.0;
  T yCm = 0.0;
  T theta = 0.0;
  T Vx = 0.0;
  T Vy = 0.0;
  T omegaZ = 0.0;
};

struct BodyReferenceState {
  T xCm = 0.0;
  T yCm = 0.0;
  T theta = 0.0;
};

struct BodyInertia {
  T bodyCenterlineLength = 0.0;
  T bodyLength = 0.0;
  T bodyWidth = 0.0;
  T bodyArea = 0.0;
  T rhoBody = 1.0;
  T mass = 0.0;
  T Ibody = 0.0;
  T semiA = 0.0;
  T semiB = 0.0;
  T mAddedSurgeTheory = 0.0;
  T mAddedHeaveTheory = 0.0;
  T IaddedTheory = 0.0;
  T mAddedSurge = 0.0;
  T mAddedHeave = 0.0;
  T Iadded = 0.0;
  T Msurge = 0.0;
  T Mheave = 0.0;
  T Itotal = 0.0;
};

struct ClampAudit {
  std::uint64_t initialPlacementClampCount = 0;
  std::uint64_t runtimeDomainClampCount = 0;
  std::uint64_t speedClampCount = 0;
  std::uint64_t omegaClampCount = 0;
  bool initialPlacementClamped = false;
  bool runtimeDomainClampHit = false;
  bool speedClampHit = false;
  bool omegaClampHit = false;
};

BodyInertia computeBodyInertia(const EelParams& p, T rhoBody = 1.0);
void enforceModeDefinitionState(SimulationCase simCase,
                                const BodyReferenceState& reference,
                                BodyState& state);
void applyRigidBodyForceUpdate(SimulationCase simCase,
                               T tKinematics,
                               T restTime,
                               T fxS,
                               T fyS,
                               T tzS,
                               const BodyInertia& inertia,
                               BodyState& state);
void advanceRigidBodyPose(SimulationCase simCase,
                          const BodyReferenceState& reference,
                          BodyState& state);
