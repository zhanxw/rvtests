#include <stdio.h>
#include <stdlib.h>

#include "base/SimpleMatrix.h"
#include "libVcf/PlinkOutputFile.h"

int main(int argc, char** argv) {
  double genotype[] = {0, 2, -9, 2, 2, 2, 2,  -9, 1,
                       2, 2, 2,  2, 1, 1, -9, -9, 0};

  SimpleMatrix m;  // marker by sample matrix
  m.resize(3, 6);
  int offset = 0;
  for (int i = 0; i < m.nrow(); ++i) {
    for (int j = 0; j < m.ncol(); ++j) {
      m[i][j] = genotype[offset];
      ++offset;
    }
  }

  std::vector<std::string> iid;
  char buf[128];
  for (int i = 0; i < m.ncol(); ++i) {
    sprintf(buf, "%d", i + 1);
    iid.push_back(buf);
  }

  std::vector<double> pheno;
  pheno.push_back(-9);
  pheno.push_back(-9);
  pheno.push_back(2);
  pheno.push_back(-9);
  pheno.push_back(2);
  pheno.push_back(2);

  PlinkOutputFile pout("testPlinkOutputFile.output");
  pout.writeFAM(iid, iid, pheno);
  pout.writeBIM("1", "snp1", 0, 1, "G", "A");
  pout.writeBIM("1", "snp2", 0, 2, "1", "2");
  pout.writeBIM("1", "snp3", 0, 3, "A", "C");
  pout.writeBED(&m, m.ncol(), m.nrow());

  fprintf(stderr, "Done\n");

  return 0;
}
