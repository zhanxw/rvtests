#include "SimpleString.h"

#include <stdlib.h>

#include <cassert>
#include <string>
#include <vector>

using std::vector;
using std::string;

char a[] = "abc\ndef\nHAH!!DFDSLIJFDSO\nadfa";
char buffer[1024];

int main(int argc, char *argv[]) {
  const std::string aa = a;
  {
    SimpleString s(4);
    for (size_t i = 0; i != aa.size(); ++i) {
      s.append(aa[i]);
    }
    assert(aa == s.data());
  }
  {
    SimpleString s(1024);
    for (size_t i = 0; i != aa.size(); ++i) {
      s.append(aa[i]);
    }
    assert(aa == s.data());
  }
}
