#include "PlinkInputFile.h"
#include "SimpleMatrix.h"

void printMatrix(SimpleMatrix& m, PlinkInputFile& f1) {
  for (int i = 0; i < m.nrow(); i++) {
    for (int j = 0; j < m.ncol(); j++) {
      int geno = (int)m[i][j];
      char ref = f1.ref[j][0];
      char alt = f1.alt[j][0];
      switch (geno) {
        case 0:
          printf("%c %c\t", ref, ref);
          break;
        case 1:
          printf("%c %c\t", ref, alt);
          break;
        case 2:
          printf("%c %c\t", alt, alt);
          break;
        default:
          if (geno < 0) printf(". .\t");
          break;
          printf("ERROR reading!\n");
      }
    }
    printf("\n");
  }
}

int main(int argc, char* argv[]) {
  // check loading all genotypes to memory
  PlinkInputFile f1("testPlinkInputFile");
  std::vector<double> maf;
  f1.calculateMAF(&maf);

  for (int i = 0; i < maf.size(); i++) {
    printf("af[%d] = %g\n", i, maf[i]);
  }

  std::vector<double> imiss;
  std::vector<double> lmiss;
  f1.calculateMissing(&imiss, &lmiss);
  for (int i = 0; i < imiss.size(); i++) {
    printf("imiss[%d] = %g\n", i, imiss[i]);
  }
  for (int i = 0; i < lmiss.size(); i++) {
    printf("lmiss[%d] = %g\n", i, lmiss[i]);
  }

  return 0;
}
