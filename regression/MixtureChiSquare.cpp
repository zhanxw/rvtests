#include "MixtureChiSquare.h"

#include <stdlib.h>
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

  // dumpLambda();
  // fprintf(stderr, "Q = %g\n", Q);
  // fprintf(stderr, "fault = %d\n", fault);
  return pValue;
};
void MixtureChiSquare::dumpLambda() const {
  for (int i = 0; i < lambda_size; ++i) {
    fprintf(stderr, "lambda[%d] = %g\n", i, lambda[i]);
  }
}
