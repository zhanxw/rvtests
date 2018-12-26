#ifndef VERSIONCHECKER_H
#define VERSIONCHECKER_H

#include <string>
#include <vector>

class VersionChecker {
 public:
  VersionChecker();
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
  /**
   * Print line 2 and onward that were retrieved earlier,
   * as usually the line 1 is software version number.
   */
  void printRemoteContent(int startLine = 1) const;

 private:
  std::string localVersion;
  std::string remoteVersion;
  std::vector<std::string> remoteInformation;
  bool quiet;
};

#endif /* VERSIONCHECKER_H */
