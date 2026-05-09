#include "physics/rigid_body.hpp"

#include <cmath>

BodyInertia computeBodyInertia(const EelParams& p, T rhoBody)
{
  BodyInertia out;
  out.bodyCenterlineLength = p.centerlineLengthLU();
  out.bodyLength = p.totalGeometricLengthLU();
  out.bodyWidth = p.bodyWidthLU();
  out.bodyArea = p.capsuleAreaLU();
  out.rhoBody = rhoBody;
  out.mass = out.bodyArea * out.rhoBody;
  out.Ibody = out.mass / 12.0 *
              (out.bodyLength * out.bodyLength + out.bodyWidth * out.bodyWidth);

  out.semiA = 0.5 * out.bodyLength;
  out.semiB = 0.5 * out.bodyWidth;
  out.mAddedSurgeTheory = out.rhoBody * M_PI * out.semiB * out.semiB;
  out.mAddedHeaveTheory = out.rhoBody * M_PI * out.semiA * out.semiA;
  const T diffAb2 = out.semiA * out.semiA - out.semiB * out.semiB;
  out.IaddedTheory = out.rhoBody * M_PI / 8.0 * diffAb2 * diffAb2;

  out.mAddedSurge = p.addedMassFrac * out.mAddedSurgeTheory;
  out.mAddedHeave = p.addedMassFrac * out.mAddedHeaveTheory;
  out.Iadded = p.addedMassFrac * out.IaddedTheory;

  out.Msurge = out.mass + out.mAddedSurge;
  out.Mheave = out.mass + out.mAddedHeave;
  out.Itotal = out.Ibody + out.Iadded;
  return out;
}

// Surge-only rigid-body model: only Vx is a dynamic state, heave/yaw are
// pinned to the reference pose by definition.
void enforceModeDefinitionState(const BodyReferenceState& reference,
                                BodyState& state)
{
  state.Vy = T(0);
  state.omegaZ = T(0);
  state.yCm = reference.yCm;
  state.theta = reference.theta;
}

void applyRigidBodyForceUpdate(T tKinematics,
                               T restTime,
                               T fxS,
                               const BodyInertia& inertia,
                               BodyState& state)
{
  if (tKinematics < restTime) {
    return;
  }
  state.Vx += fxS / inertia.Msurge;
}

void advanceRigidBodyPose(const BodyReferenceState& reference,
                          BodyState& state)
{
  state.xCm += state.Vx;
  state.yCm = reference.yCm;
  state.theta = reference.theta;
}
