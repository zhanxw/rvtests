#include <stdio.h>

#include "BCFReader.h"
#include "VCFExtractor.h"

int main(int argc, char* argv[]) {
  {
    std::string fn = "all.anno.filtered.extract.bcf.gz";
    std::string r = "1:196621007-196716634";
    VCFExtractor vin(fn);
    vin.setRangeList(r);
    printf("Extract two sites\n");
    vin.setAnnoType("Nonsynonymous");
    while (vin.readRecord()) {
      VCFRecord& r = vin.getVCFRecord();
      VCFPeople& people = r.getPeople();
      VCFIndividual* indv;

      printf("%s:%d\n", r.getChrom(), r.getPos());
    }
  }

#if 0
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
#endif

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
    assert(n == 14);  // without header lines
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
    assert(n == 2);
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
    std::string h = tr.getHeader().c_str();
    int count = 0;
    for (size_t i = 0; i < h.size(); ++i) {
      if (h[i] == '\n') count++;
    }
    fprintf(stdout, "header has %d lines.\n", count);
    assert(count == 60);
  }
  return 0;
}
