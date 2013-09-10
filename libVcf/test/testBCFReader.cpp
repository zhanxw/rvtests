#include "BCFReader.h"

int main(int argc, char *argv[])
{
  const char* fn = "test.bcf.gz";
  {
    BCFReader tr(fn);
    int n = 0;
    std::string line;
    while (tr.readLine(&line)) {
      // printf("line %d = [ %s ]\n", n, line.c_str());
      ++n;
    }
    fprintf(stdout, "Read %d lines\n", n);
    assert ( n == 14) ; // without header lines
  }
  {
    BCFReader tr(fn);
    tr.addRange("1:196341181-196341254");
    int n = 0;
    std::string line;
    while (tr.readLine(&line)) {
      // printf("line %d = [ %s ]\n", n, line.c_str());
      ++n;
    }
    fprintf(stdout, "Read %d lines\n", n);
    assert ( n == 2) ;
  }
  {
    BCFReader tr(fn);
    tr.addRange("1:196341857-196341857");
    tr.addRange("1:196341181-196341254");
    int n = 0;
    std::string line;
    while (tr.readLine(&line)) {
      // printf("line %d = [ %s ]\n", n, line.c_str());
      ++n;
    }
    fprintf(stdout, "Read %d lines\n", n);
    // assert ( n == 3) ;
  }

  {
    BCFReader tr(fn);
    tr.addRange("2:196341857-196341857");
    tr.addRange("2:196341181-196341254");
    int n = 0;
    std::string line;
    while (tr.readLine(&line)) {
      printf("line %d = [ %s ]\n", n, line.c_str());
      ++n;
    }
    fprintf(stdout, "Read %d lines\n", n);
    // assert ( n == 0) ;
  }

  return 0;
}
