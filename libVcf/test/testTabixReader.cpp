#include "TabixReader.h"

int main(int argc, char *argv[])
{
  const char* fn = "test.vcf.gz";
  {
    TabixReader tr(fn);
    int n = 0;
    std::string line;
    while (tr.readLine(&line)) {
      // printf("line %d = [ %s ]\n", n, line.c_str());
      ++n;
    }
    fprintf(stdout, "Read %d lines\n", n);
    assert ( n == 14) ;
  }

  {
    TabixReader tr(fn);
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
    TabixReader tr(fn);
    tr.addRange("1:196341857-196341857");
    tr.addRange("1:196341181-196341254");
    int n = 0;
    std::string line;
    while (tr.readLine(&line)) {
      // printf("line %d = [ %s ]\n", n, line.c_str());
      ++n;
    }
    fprintf(stdout, "Read %d lines\n", n);
    assert ( n == 3) ;
  }

  {
    TabixReader tr(fn);
    tr.addRange("2:196341857-196341857");
    tr.addRange("2:196341181-196341254");
    int n = 0;
    std::string line;
    while (tr.readLine(&line)) {
      printf("line %d = [ %s ]\n", n, line.c_str());
      ++n;
    }
    fprintf(stdout, "Read %d lines\n", n);
    assert ( n == 0) ;
  }

  {
    TabixReader tr(fn);
      std::string h= tr.getHeader().c_str();
    int count = 0;
    for (size_t i = 0; i < h.size(); ++i) {
      if (h[i] == '\n')
        count++;
    }
    fprintf(stdout, "header has %d lines.\n", count);
    assert(count == 36);
  }
  return 0;
}
