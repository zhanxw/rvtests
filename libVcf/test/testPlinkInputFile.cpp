#include "PlinkInputFile.h"
#include "SimpleMatrix.h"

void printMatrix(SimpleMatrix& m, PlinkInputFile& f1) {
  for (int i = 0; i < m.nrow(); i++ ) {
    for (int j = 0; j < m.ncol(); j++) {
      int geno = (int) m[i][j];
      char ref = f1.ref[j];
      char alt = f1.alt[j];
      switch(geno) {
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
          if ( geno < 0)
            printf(". .\t");
          break;
          printf("ERROR reading!\n");
      }
    }
    printf("\n");
  }
}

int main(int argc, char *argv[])
{
  // check loading all genotypes to memory
  PlinkInputFile f1 ("testPlinkInputFile");
  SimpleMatrix m;
  f1.readIntoMatrix(&m);
  for (int i = 0; i < m.nrow(); i++ ) {
    for (int j = 0; j < m.ncol(); j++) {
      int geno = (int) m[i][j];
      printf("%d ", geno);
    }
    printf("\n");
  }
  printMatrix(m, f1);

  // loading some people and marker
  printf("----------------------------------------------------------------------\n");
  m.clear();
  f1.readIntoMatrix(&m, NULL, NULL);
  std::vector<std::string> personName, markerName;
  f1.readIntoMatrix(&m, &personName, &markerName);      
  printMatrix(m, f1);


  // loading some people and marker
  printf("----------------------------------------------------------------------\n");
  m.clear();
  std::vector<std::string> people;
  std::vector<std::string> marker;
  people.push_back("1");
  f1.readIntoMatrix(&m, &people, &marker);
  printMatrix(m, f1);


  // loading some people and marker
  printf("----------------------------------------------------------------------\n");
  m.clear();
  people.push_back("6");
  marker.push_back("snp1");
  f1.readIntoMatrix(&m, &people, &marker);
  printMatrix(m, f1);

  // loading some people and marker
  printf("----------------------------------------------------------------------\n");
  m.clear();
  people.push_back("2");
  marker.push_back("snp2");
  f1.readIntoMatrix(&m, &people, &marker);
  printMatrix(m, f1);

  return 0;
}
