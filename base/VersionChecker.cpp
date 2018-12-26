#include "VersionChecker.h"

#include "Http.h"
#include "TypeConversion.h"

VersionChecker::VersionChecker() {
  this->quiet = (std::getenv("VERSIONCHECKER_DEBUG") == NULL);
}

int VersionChecker::retrieveRemoteVersion(const std::string& urlToVersion) {
  Http http(urlToVersion, 0.5);
  HttpResponse response;
  if (this->quiet) {
    http.enableQuiet();
  } else {
    http.disableQuiet();
  }
  if (!this->quiet) {
    fprintf(stderr, "retrieve remote version from [ %s ]\n",
            urlToVersion.c_str());
  }
  // 2: the remote version file has two lines
  int ret = http.read(&response, 2);
  if (ret < 0) {
    if (!this->quiet) {
      fprintf(
          stderr,
          "Failed to retrieve remote version, probably due to wrong http\n");
    }
    return -1;
  }
  this->remoteInformation = response.getBody();
  if (this->remoteInformation.size() < 1) {
    if (!this->quiet) {
      fprintf(stderr,
              "Failed to retrieve remote version due to empty version "
              "information\n");
    }
    return -1;
  }
  this->remoteVersion = this->remoteInformation[0];
  if (!this->quiet) {
    fprintf(stderr, "Retrieved remote version [ %s ]\n",
            this->remoteVersion.c_str());
  }
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
