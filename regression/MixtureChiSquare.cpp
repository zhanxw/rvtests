#include "MixtureChiSquare.h"
#include "qfc.c"

double MixtureChiSquare::getPvalue(double Q) {
  int fault;
  double trace[7];

  double pValue = 1.0 - qf(lambda, noncen, df, lambda_size, sigma,
                           Q, lim, acc, trace, &fault);
  if(pValue>1.0)
    pValue = 1.0; //this occurs when eigen values are very large
  if (fault) {
    pValue = -1.0; //
  }
  return pValue;
};
