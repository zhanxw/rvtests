#ifndef VERSIONCHECKER_H
#define VERSIONCHECKER_H

#include <string>
#include <vector>

class VersionChecker {
 public:
  /**
   * @return 0 only if succeed 
   */
  int retrieveRemoteVersion(const std::string& urlToRemoteVersion);
  int setLocalVersion(const std::string& localVersion);
  bool isRemoteVersionNewer() const;
  /**
   * Check if remote version is newer than the @param currentVersion
   * @return 1 if remote is newer; or 0 if not newer; or -1 if error occurs
   */
  const std::vector<std::string>& getRemoteContent() const {
    return this->remoteInformation;
  }
  void printRemoteContent() const;
 private:
  std::string localVersion;
  std::string remoteVersion;
  std::vector<std::string> remoteInformation;
};

#endif /* VERSIONCHECKER_H */
