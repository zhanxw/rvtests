#ifndef _KINSHIP_H_
#define _KINSHIP_H_

#include "Pedigree.h"
#include "SimpleMatrix.h"

namespace zhanxw{

  class Kinship {
 public:
    int constructFromPedigree(const zhanxw::Pedigree& ped);
    void dumpKinshipFromPedigree(const zhanxw::Pedigree& ped) {
      int n = ped.getPeopleNumber();
      for (int i = 0; i < n; ++i) {
        const zhanxw::Person& p = ped.getPeople()[i];
        printf("%s\t%s", ped.getFamilyName(p.getFamily()), ped.getPersonName(i));
        for (int j = 0; j < n; ++j) {
          printf("\t%g", mat[i][j]);
        }
        printf("\n");
      }
    };
    int constructFromGenotype(const SimpleMatrix& m);
    
 private:
    SimpleMatrix mat;
  };
} // namespace

#endif /* _KINSHIP_H_ */
