//#include "Pedigree.h"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string>
#include "base/Pedigree.h"
#include "base/Kinship.h"

void usage(int argc, char** argv) {
  fprintf(stderr, "%s [a|x] pedFile listFile\n", argv[0]);
  exit(1);
}

bool isSameFamily(const zhanxw::Pedigree&ped, int i, int j) {
  return (ped.getPeople()[i].getFamily() == ped.getPeople()[j].getFamily());
}

/**
 * @param fn has two columns, which are the fam id an person id to include
 */
void outputKinship(const char* fn, const zhanxw::Pedigree& ped, const SimpleMatrix& m) {
  std::map< std::string, double> result;
  int n = m.nrow();
  std::string key;

  std::set<std::string> s;
  LineReader lr(fn);
  std::vector<std::string> fd;
  while(lr.readLineBySep(&fd, " \t")) {
    key.clear();
    // key += fd[0];
    // key += "->";
    key += fd[0];
    s.insert(key);
  }
  
  for (int i = 0; i < n; ++i) {
    for (int j = i; j < n; ++j) {
      key.clear();
      if (s.count(ped.getPersonName(i)) &&
          s.count(ped.getPersonName(j))
          ) {
        int f = ped.getPeople()[i].getFamily();
        std::string fam = ped.getFamilyName(f);
        if (!isSameFamily(ped, i, j)) continue;
        printf("%s %s %s %.6f\n",
               fam.c_str(),
               ped.getPersonName(i),
               ped.getPersonName(j),
               m[i][j]);
      }
    }
  }

  
  //   double val = -1;
  //   if (result.count(key)) {
  //     val = result[key];
  //   }
  //   fprintf(stdout, "%s %s %s %g\n", fd[0].c_str(),
  //           fd[1].c_str(), fd[2].c_str(), val);
  // }
}

int main(int argc, char *argv[])
{
  if (argc != 4) {
    usage(argc, argv);
  }

  zhanxw::Pedigree ped;
  int ret;
  // "testKinship.ped"
  ret = loadPedigree(argv[2], &ped);
  if (!ret) {
    // printf("Loaded OK\n");
    // printf("Loaded %zu family and %zu people.\n", ped.getFamilyNumber(), ped.getPeopleNumber());
  } else {
    fprintf(stderr, "Load failed!\n");
  }

  if (argv[1][0] == 'a') { // check autosomal
    zhanxw::Kinship kin;
    kin.constructFromPedigree(ped);
    outputKinship(argv[3], ped, kin.getKinship());
  } else if (argv[1][0] == 'x') { // check X chrom
    zhanxw::KinshipForX kin;
    kin.constructFromPedigree(ped);
    outputKinship(argv[3], ped, kin.getKinship());
  }
  
  return 0;
}

