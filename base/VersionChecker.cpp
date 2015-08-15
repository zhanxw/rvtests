#include "VersionChecker.h"
#include "Http.h"
#include "TypeConversion.h"


VersionChecker::VersionChecker(const std::string& urlToVersion) {
  Http http(urlToVersion);
  http.read(&this->remoteInformation);
}

int VersionChecker::hasNewVersionThan(const std::string& currentVersion) {
  int ver;
  if (str2int(currentVersion, &ver)) {
    return 0;
  }
  return this->hasNewVersionThan(ver);
}

int VersionChecker::hasNewVersionThan(int currentVersion) {
  if (!this->remoteInformation.size()) {
    return 0;
  }
  int remoteVersion;
  if (!str2int(remoteInformation[0], &remoteVersion)) return 0;
  if (remoteVersion > currentVersion) return 1;
  return 0;
}

void VersionChecker::printNewVersion() const {
  if (!remoteInformation.size()) {
    return;
  }
  for (size_t i = 1; i != remoteInformation.size(); ++i) {
    fprintf(stderr, "%s\n", remoteInformation[i].c_str());
  }
}
