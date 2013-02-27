//#include "Pedigree.h"
#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <string>
#include "base/Pedigree.h"

int main(int argc, char *argv[])
{
  zhanxw::Pedigree ped;
  int ret;

  ret = loadPedigree(argv[1], &ped);
  if (!ret) {
    printf("Loaded OK\n");
    printf("Loaded %zu family and %zu people.\n", ped.getFamilyNumber(), ped.getPeopleNumber());
  } else {
    fprintf(stderr, "Load failed!\n");
  }

  dumpPedigree(ped);
  dumpPedFile(ped);
  
  return 0;
}
