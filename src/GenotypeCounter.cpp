#include "GenotypeCounter.h"

#include <stdio.h>
#include <stdlib.h>
#include "snp_hwe.c"

double GenotypeCounter::getHWE() const {
  double hweP = 0.0;
  if (nHomRef + nHet + nHomAlt == 0 ||
      (nHet < 0 || nHomRef < 0 || nHomAlt < 0)) {
    hweP = 0.0;
  } else {
    hweP = SNPHWE(nHet, nHomRef, nHomAlt);
  }
  return hweP;
}
