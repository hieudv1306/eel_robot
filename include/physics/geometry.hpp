#pragma once

#include "core/params.hpp"
#include "core/types.hpp"

#include <array>
#include <vector>

void buildCapsuleGeometry(
    const EelParams& p,
    T t, T xCm, T yCm, T theta,
    T dtLbm, T ampRamp,
    std::vector<T>& globX, std::vector<T>& globY,
    std::vector<T>& vxDefGlob, std::vector<T>& vyDefGlob,
    std::vector<T>& allDs);

std::array<T, 2> bodyForwardAxis(T theta);
std::array<T, 2> bodyLateralAxis(T theta);
