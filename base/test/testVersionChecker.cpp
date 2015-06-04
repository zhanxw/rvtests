#include "VersionChecker.h"

int main(int argc, char *argv[])
{
  {
    fprintf(stderr, "Has new version\n");
    VersionChecker ver("http://zhanxw.com/rvtests/version");
    if (ver.hasNewVersionThan(-1)) {
      ver.printNewVersion();
    }
  }
  {
    fprintf(stderr, "No new version\n");
    VersionChecker ver("http://zhanxw.com/rvtests/version");
    if (ver.hasNewVersionThan(99999999)) {
      ver.printNewVersion();
    }
  }
  {
    fprintf(stderr, "Wrong address\n")
        ;
    VersionChecker ver("http://nonexist.zhanxw/rvtests/version");
    if (ver.hasNewVersionThan(99999999)) {
      ver.printNewVersion();
    }
  }
  return 0;
}
