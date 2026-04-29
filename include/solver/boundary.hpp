#pragma once

#include "core/params.hpp"
#include "core/types.hpp"
#include "solver/openlb_setup.hpp"

inline void resetEulerForceField(
    const EelParams& p,
    olb::SuperLattice<T, DESCRIPTOR>& sLattice,
    bool fixedInflow)
{
  const int nx = p.nx;
  const int ny = p.ny;
  const int spW = p.spongeWidth;
  const T spS = p.spongeStrength;

  sLattice.setProcessingContext(olb::ProcessingContext::Evaluation);
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      T force[2] = {0.0, 0.0};
      auto cell = sLattice.get(0, i, j);
      cell.template setField<olb::descriptors::FORCE>({0.0, 0.0});
      T sigma = 0.0;
      if (!fixedInflow && spW > 0 && spS > 0.0 && i < spW) {
        T xi = T(spW - i) / T(spW);
        sigma = spS * xi * xi;
      }
      else if (spW > 0 && spS > 0.0 && i >= nx - spW) {
        T xi = T(i - (nx - spW - 1)) / T(spW);
        sigma = spS * xi * xi;
      }
      if (sigma > 0.0) {
        T uCell[2] = {0.0, 0.0};
        cell.computeU(uCell);
        const T targetUx = fixedInflow ? p.inflowVelocity : T(0);
        force[0] = -sigma * (uCell[0] - targetUx);
        force[1] = -sigma * uCell[1];
      }
      cell.template setField<olb::descriptors::FORCE>({force[0], force[1]});
    }
  }
  sLattice.template setProcessingContext<
    olb::Array<olb::descriptors::FORCE>>(olb::ProcessingContext::Simulation);
}
