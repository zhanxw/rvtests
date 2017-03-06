#include "IO.h"
#include "Utils.h"

#include <stdlib.h>

#include <cassert>
#include <string>
#include <vector>

using std::vector;
using std::string;

char a[] = "abc\ndef\nHAH!!DFDSLIJFDSO\nadfa";
char buffer[1024];

int main(int argc, char *argv[]) {
  vector<string> t;
  int ret = stringTokenize(a, '\n', &t);
  assert(ret == 4);

  {
    char fn[] = "http://zhanxw.com/rvtests/version";
    LineReader lr(fn);
    std::string ln;
    int i = 0;
    while (lr.readLine(&ln)) {
      fprintf(stderr, "line %d: %s\n", i + 1, ln.c_str());
      i++;
    }
  }
  {
    char fn[] = "http://zhanxw.com/rvtests/version.gz";
    LineReader lr(fn);
    std::string ln;
    int i = 0;
    while (lr.readLine(&ln)) {
      fprintf(stderr, "line %d: %s\n", i + 1, ln.c_str());
      i++;
    }
  }

  return 0;
}
