#pragma once

#include "core/params.hpp"
#include "core/types.hpp"
#include "physics/soft_backbone.hpp"

#include <array>
#include <vector>

void buildCapsuleGeometry(
    const EelParams& p,
    T t, T xCm, T yCm, T theta,
    T dtLbm, T ampRamp,
    std::vector<T>& globX, std::vector<T>& globY,
    std::vector<T>& vxDefGlob, std::vector<T>& vyDefGlob,
    std::vector<T>& allDs);

void buildCapsuleGeometryFromBackboneState(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    T xCm, T yCm, T theta,
    std::vector<T>& globX, std::vector<T>& globY,
    std::vector<T>& allDs);

void buildCapsuleGeometryFromBackboneMotion(
    const EelParams& p,
    const SoftBackboneConfig& config,
    const SoftBackboneState& state,
    const SoftBackboneState& statePlus,
    const SoftBackboneState& stateMinus,
    T xCm, T yCm, T theta,
    std::vector<T>& globX, std::vector<T>& globY,
    std::vector<T>& vxDefGlob, std::vector<T>& vyDefGlob,
    std::vector<T>& allDs);

std::array<T, 2> bodyForwardAxis(T theta);
std::array<T, 2> bodyLateralAxis(T theta);
