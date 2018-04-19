#include <stdio.h>
#include <stdlib.h>

#include "KGGInputFile.h"
#include "base/Utils.h"

void printVCFMeta() {
  // GP is between 0 and 1 in VCF v4.3, but phred-scaled value in VCF v4.2
  printf("##fileformat=VCFv4.3\n");
  printf(
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype "
      "calls\">\n");
}

void printVCFHeader(const std::vector<std::string>& sm) {
  printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  for (size_t i = 0; i != sm.size(); ++i) {
    printf("\t%s", sm[i].c_str());
  }
  printf("\n");
}

void printGeno(int i, const std::vector<int>& g) {
  printf("%d => ", i);
  for (int i = 0; i < g.size(); ++i) {
    printf(" %d", g[i]);
  }
  printf("\n");
}

int main(int argc, char* argv[]) {
  KGGInputFile read("kggseq1", "gz");

  const int N = read.getNumSample();
  const int M = read.getNumMarker();
  const std::vector<std::string> sm = read.getSampleName();
  printVCFMeta();
  printVCFHeader(sm);

  int a1, a2;
  for (int i = 0; i < M; ++i) {
    if (!read.readRecord()) {
      continue;
    }
    printf("%s", read.getChrom()[i].c_str());
    printf("\t%d", read.getPosition()[i]);
    printf("\t%s", read.getMarkerName()[i].c_str());
    printf("\t%s", read.getRef()[i].c_str());
    printf("\t%s", stringJoin(read.getAlt()[i], ',').c_str());
    printf("\t.\t.\t.\tGT");
    for (int j = 0; j < N; ++j) {
      read.getAllele(j, &a1, &a2);
      printf("\t%d/%d", a1, a2);
    }
    printf("\n");
  }  // loop marker

  return 0;
}
