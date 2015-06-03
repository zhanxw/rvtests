#ifndef VERSIONCHECKER_H
#define VERSIONCHECKER_H

#include <string>
#include <vector>

class VersionChecker {
 public:
  VersionChecker(const std::string& urlToVersion);
  int hasNewVersionThan(const std::string& currentVersion);  
  int hasNewVersionThan(int currentVersion);
  void printNewVersion() const;
 protected:
  std::vector<std::string> remoteInformation;
};

#endif /* VERSIONCHECKER_H */
