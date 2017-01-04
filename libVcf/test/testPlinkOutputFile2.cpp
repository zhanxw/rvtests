#include <stdio.h>
#include <stdlib.h>

#include "base/SimpleMatrix.h"
#include "libVcf/PlinkOutputFile.h"

std::vector<int> seq(int n) {
  std::vector<int> ret(n);
  for (int i = 0; i < n; ++i) {
    ret[i] = i;
  }
  return ret;
}

int main(int argc, char** argv) {
  PlinkOutputFile pout("testPlinkOutputFile.output");

  std::vector<int> sampleIdx = seq(6);
  std::vector<int> markerIdx = seq(3);

  pout.extract("testPlinkInputFile", sampleIdx, markerIdx);

  fprintf(stderr, "Done\n");

  return 0;
}
