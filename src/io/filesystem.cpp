#include "io/filesystem.hpp"

#include <cerrno>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

bool eelDirectoryExists(const std::string& path);

bool eelEnsureDirectory(const std::string& path)
{
  if (path.empty()) {
    return false;
  }
  if (::mkdir(path.c_str(), 0755) == 0) {
    return true;
  }
  if (errno == EEXIST) {
    return eelDirectoryExists(path);
  }
  return false;
}

bool eelDirectoryExists(const std::string& path)
{
  struct stat st;
  return ::stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

std::string trimTrailingSlash(const std::string& path)
{
  if (path.size() <= 1) {
    return path;
  }
  size_t end = path.size();
  while (end > 1 && path[end - 1] == '/') {
    --end;
  }
  return path.substr(0, end);
}

bool eelEnsureDirectoryTree(const std::string& rawPath)
{
  if (rawPath.empty()) {
    return false;
  }

  const std::string path = trimTrailingSlash(rawPath);
  std::string partial;
  size_t pos = 0;
  if (!path.empty() && path[0] == '/') {
    partial = "/";
    pos = 1;
  }

  while (pos <= path.size()) {
    const size_t next = path.find('/', pos);
    const std::string part = path.substr(pos, next - pos);
    if (!part.empty()) {
      if (!partial.empty() && partial.back() != '/') {
        partial += "/";
      }
      partial += part;
      if (!eelEnsureDirectory(partial)) {
        return false;
      }
    }
    if (next == std::string::npos) {
      break;
    }
    pos = next + 1;
  }
  return true;
}

std::string sanitizePathToken(const std::string& raw)
{
  std::string out;
  out.reserve(raw.size());
  bool lastWasUnderscore = false;
  for (unsigned char ch : raw) {
    char c = static_cast<char>(ch);
    if (std::isalnum(ch) || c == '-') {
      out.push_back(c);
      lastWasUnderscore = false;
    } else if (c == '.') {
      out.push_back('p');
      lastWasUnderscore = false;
    } else if (!lastWasUnderscore) {
      out.push_back('_');
      lastWasUnderscore = true;
    }
  }
  while (!out.empty() && out.back() == '_') {
    out.pop_back();
  }
  return out.empty() ? "run" : out;
}

std::string formatValueForRunId(T value)
{
  std::ostringstream os;
  os << std::setprecision(8) << std::defaultfloat << value;
  std::string s = os.str();
  for (char& c : s) {
    if (c == '-') {
      c = 'm';
    } else if (c == '+') {
      c = 'p';
    }
  }
  return sanitizePathToken(s);
}

std::string appendSuffixBeforeExtension(const std::string& path,
                                        const std::string& suffix)
{
  const size_t slashPos = path.find_last_of("/\\");
  const size_t dotPos = path.find_last_of('.');
  if (dotPos == std::string::npos ||
      (slashPos != std::string::npos && dotPos < slashPos)) {
    return path + suffix;
  }
  return path.substr(0, dotPos) + suffix + path.substr(dotPos);
}

std::string makeBaseRunId(const EelParams& p, SimulationCase simCase,
                          const std::string& runTag)
{
  std::string runId = std::string("case_") + simulationCaseName(simCase)
                    + "_AR" + formatValueForRunId(p.aspectRatio)
                    + "_tau" + formatValueForRunId(p.tau)
                    + "_kap" + formatValueForRunId(p.kappa)
                    + "_sub" + std::to_string(p.substeps);
  runId += "_N" + std::to_string(p.nx) + "x" + std::to_string(p.ny);
  if (!runTag.empty()) {
    runId += "_" + sanitizePathToken(runTag);
  }
  return runId;
}

std::string makeUniqueRunId(const std::string& baseRunId,
                            const std::string& runsRootDir)
{
  std::string candidate = baseRunId;
  for (int repeat = 2; eelDirectoryExists(runsRootDir + candidate); ++repeat) {
    candidate = baseRunId + "_r" + std::to_string(repeat);
  }
  return candidate;
}

