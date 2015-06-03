#include "VersionChecker.h"
#include "IO.h"
#include "TypeConversion.h"

VersionChecker::VersionChecker(const std::string& urlToVersion) {
#ifdef _USE_KNETFILE
  LineReader lr(urlToVersion);
  std::string line;
  while (lr.readLine(&line)) {
    this->remoteInformation.push_back(line);
  }
#endif
}

int VersionChecker::hasNewVersionThan(const std::string& currentVersion) {
  int ver;
  if (str2int(currentVersion, &ver)) {
    return this->hasNewVersionThan(ver);
  }
  return 0;
}

int VersionChecker::hasNewVersionThan(int currentVersion) {
  int remoteVersion;
  if (!str2int(remoteInformation[0], &remoteVersion)) return 0;
  if (remoteVersion > currentVersion) return 1;
  return 0;
}

void VersionChecker::printNewVersion() const {
  for (size_t i = 1; i != remoteInformation.size(); ++i) {
    fprintf(stderr, "%s\n", remoteInformation[i].c_str());
  }
}
