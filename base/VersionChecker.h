#ifndef VERSIONCHECKER_H
#define VERSIONCHECKER_H

#include <string>
#include <vector>

class VersionChecker {
 public:
  VersionChecker(const std::string& urlToVersion);
  /**
   * Check if remote version is newer than the @param currentVersion
   * @return 1 if remote is newer; or 0 if not newer or if error occurs
   */
  int hasNewVersionThan(const std::string& currentVersion);  
  int hasNewVersionThan(int currentVersion);
  void printNewVersion() const;
 protected:
  std::vector<std::string> remoteInformation;
};

#endif /* VERSIONCHECKER_H */
