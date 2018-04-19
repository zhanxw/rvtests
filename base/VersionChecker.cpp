#include "VersionChecker.h"

#include "Http.h"
#include "TypeConversion.h"

int VersionChecker::retrieveRemoteVersion(const std::string& urlToVersion) {
  Http http(urlToVersion, 0.5);
  http.enableQuiet();

  if (http.read(&this->remoteInformation) < 0) {
    return -1;
  }
  if (this->remoteInformation.size() < 1) {
    return -1;
  }
  this->remoteVersion = this->remoteInformation[0];
  return 0;
}

int VersionChecker::setLocalVersion(const std::string& currentVersion) {
  this->localVersion = currentVersion;
  return 0;
}

bool VersionChecker::isRemoteVersionNewer() const {
  int local = 0;
  int remote = 0;
  if (!str2int(this->localVersion, &local)) return false;
  if (!str2int(this->remoteVersion, &remote)) return false;
  if (remote > local) return true;
  return false;
}

void VersionChecker::printRemoteContent(int startLine) const {
  if ((int)remoteInformation.size() <= startLine) return;
  for (size_t i = startLine; i != remoteInformation.size(); ++i) {
    fprintf(stderr, "%s\n", remoteInformation[i].c_str());
  }
}
