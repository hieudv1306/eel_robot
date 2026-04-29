#pragma once

#include "core/enums.hpp"
#include "core/params.hpp"
#include "core/types.hpp"

#include <string>

bool eelDirectoryExists(const std::string& path);
bool eelEnsureDirectory(const std::string& path);
std::string trimTrailingSlash(const std::string& path);
bool eelEnsureDirectoryTree(const std::string& rawPath);
std::string sanitizePathToken(const std::string& raw);
std::string formatValueForRunId(T value);
std::string appendSuffixBeforeExtension(const std::string& path, const std::string& suffix);
std::string makeBaseRunId(const EelParams& p, SimulationCase simCase, const std::string& runTag);
std::string makeUniqueRunId(const std::string& baseRunId, const std::string& runsRootDir);
