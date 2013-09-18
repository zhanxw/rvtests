#include "BCFReader.h"

int main(int argc, char *argv[])
{
  if (false) {
    // this code demonstrate how to temporarily close stdout
    int fd = STDOUT_FILENO;
    int dupFd = dup(fileno(stdout));
    stdout = fdopen(dupFd, "w");
    assert(stdout);

    assert(dupFd > 0);
    printf("before stdout, fd = %d, dupFd = %d\n", fd, dupFd);
    fclose(stdout);
    printf("after stdout\n");
    return 0;
  }
  
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

  {
    BCFReader tr(fn);
    fprintf(stdout, "header = [%s]\n", tr.getHeader().c_str());
  }
  return 0;
}
