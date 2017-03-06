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
    char fn[] = "abc.txt";
    FileWriter fw(fn);
    fw.write(a);
    fw.close();

    LineReader lr(fn);
    std::string ln;
    int i = 0;
    while (lr.readLine(&ln)) {
      assert(0 == strcmp(t[i].c_str(), ln.c_str()));
      i++;
    }
  }

  {
    char fn[] = "abc.txt.gz";
    FileWriter fw(fn);
    fw.write(a);
    fw.close();

    LineReader lr(fn);
    std::string ln;
    int i = 0;
    while (lr.readLine(&ln)) {
      assert(0 == strcmp(t[i].c_str(), ln.c_str()));
      i++;
    }
  }

  {
    char fn[] = "abc.txt.bz2";
    FileWriter fw(fn);
    fw.write(a);
    fw.close();

    LineReader lr(fn);
    std::string ln;
    int i = 0;
    while (lr.readLine(&ln)) {
      assert(0 == strcmp(t[i].c_str(), ln.c_str()));
      i++;
    }
  }

  {
    char fn[] = "abc.txt.bgzip.gz";
    FileWriter fw(fn, BGZIP);
    fw.write(a);
    fw.close();

    LineReader lr(fn);
    std::string ln;
    int i = 0;
    while (lr.readLine(&ln)) {
      assert(0 == strcmp(t[i].c_str(), ln.c_str()));
      i++;
    }
  }

  {
    char fn[] = "abc.txt";
    FileWriter fw(fn);
    fw.write(a);
    fw.close();

    LineReader lr(fn);
    std::vector<std::string> fd;
    int i = 0;
    while (lr.readLineBySep(&fd, "\t")) {
      assert(0 == strcmp(t[i].c_str(), fd[0].c_str()));
      i++;
    }
  }

  return 0;
}
