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

  {
    char a[] =
        "r1c1 r1c2\n"
        "r2c1\tr2c2\n"
        "r3c1  r3c2\n"
        "r4c1\t\tr4c2\n"
        "r5c1 \t  r5c2 \n"
        "r6c1 r6c2";
    char fn[] = "abc.txt";
    FileWriter fw(fn);
    fw.write(a);
    fw.close();

    LineReader lr(fn);
    std::vector<std::string> fd;
    assert(lr.readLineBySep(&fd, " \t"));
    assert(fd.size() == 2);
    assert(fd[0] == "r1c1");
    assert(fd[1] == "r1c2");

    assert(lr.readLineBySep(&fd, " \t"));
    assert(fd.size() == 2);
    assert(fd[0] == "r2c1");
    assert(fd[1] == "r2c2");

    assert(lr.readLineBySep(&fd, " \t"));
    assert(fd.size() == 3);
    assert(fd[0] == "r3c1");
    assert(fd[2] == "r3c2");

    assert(lr.readLineBySep(&fd, " \t"));
    assert(fd.size() == 3);
    assert(fd[0] == "r4c1");
    assert(fd[2] == "r4c2");

    assert(lr.readLineBySep(&fd, " \t"));
    assert(fd.size() == 6);
    assert(fd[0] == "r5c1");
    assert(fd[4] == "r5c2");

    assert(lr.readLineBySep(&fd, " \t"));
    assert(fd.size() == 2);
    assert(fd[0] == "r6c1");
    assert(fd[1] == "r6c2");

    assert(!lr.readLineBySep(&fd, " \t"));
  }
  return 0;
}
