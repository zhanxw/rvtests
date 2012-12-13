//#include "Pedigree.h"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string>
#include "base/Pedigree.h"
#include "base/Kinship.h"

int main(int argc, char *argv[])
{
  zhanxw::Pedigree ped;
  int ret;
  // "testKinship.ped"
  ret = loadPedigree(argv[1], &ped);
  if (!ret) {
    printf("Loaded OK\n");
    printf("Loaded %zu family and %zu people.\n", ped.getFamilyNumber(), ped.getPeopleNumber());
  } else {
    fprintf(stderr, "Load failed!\n");
  }

  dumpPedigree(ped);
  // dumpPedFile(ped);
  std::vector<int> seq;
  if (ped.calculateIterationSequence(&seq)) {
    fprintf(stderr, "Cannot get sequenced!\n");
  }  else {
    for (size_t i = 0; i < seq.size(); ++i) {
      printf("%d ", seq[i]);
    }
    printf("\n");
  }
  
  zhanxw::Kinship kin;
  kin.constructFromPedigree(ped);
  kin.dumpKinshipFromPedigree(ped);
  
  return 0;
}

