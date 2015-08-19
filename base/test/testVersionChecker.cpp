#include "VersionChecker.h"
#include <cassert>

int main(int argc, char *argv[])
{
  {
    fprintf(stderr, "Has new version\n");
    VersionChecker ver;
    assert(ver.retrieveRemoteVersion("http://zhanxw.com/rvtests/version") == 0);
    assert(ver.setLocalVersion("-1") == 0);
    assert(ver.isRemoteVersionNewer());
    fprintf(stderr, "-- Remote Content BEG --\n");
    ver.printRemoteContent();
    fprintf(stderr, "-- Remote Content END --\n");
    fprintf(stderr, "\n");
  }
  {
    fprintf(stderr, "No new version\n");
    VersionChecker ver;
    assert(ver.retrieveRemoteVersion("http://zhanxw.com/rvtests/version") == 0);
    assert(ver.setLocalVersion("99999999") == 0);
    assert(!ver.isRemoteVersionNewer());
    fprintf(stderr, "\n");    
  }
  {
    fprintf(stderr, "Wrong address\n");
    VersionChecker ver;
    assert(ver.retrieveRemoteVersion("http://nonexist.zhanxw/rvtests/version"));
    assert(!ver.isRemoteVersionNewer()); 
    fprintf(stderr, "-- Remote Content BEG --\n");
    ver.printRemoteContent(0);
    fprintf(stderr, "-- Remote Content END --\n");
    fprintf(stderr, "\n");
  }
  return 0;
}
