#include "MixtureChiSquare.h"

#include <stdlib.h>
#include "cdflib.h" // for noncentral chi-sq
#include "qfc.c"    // for mixutre chi-sq

double MixtureChiSquare::getPvalue(double Q) {
  int fault;
  double trace[7];

  double pValue = 1.0 - qf(lambda, noncen, df, lambda_size, sigma, Q, lim, acc,
                           trace, &fault);
  if (pValue > 1.0)
    pValue = 1.0;  // this occurs when eigen values are very large
  if (fault) {
    pValue = -1.0;  //
  }

  // dumpLambda();
  // fprintf(stderr, "Q = %g\n", Q);
  // fprintf(stderr, "fault = %d\n", fault);
  return pValue;
}

double sum(double* d, int n, int power) {
  double r = 0.0;
  double tmp;
  for (int i = 0; i < n; ++i) {
    tmp = d[i];
    for (int j = 1; j < power; ++j) {
      tmp *= d[i];
    }
    r += tmp;
  }
  return r;
}

double MixtureChiSquare::getLiuPvalue(double Q) {
  const double c1 = sum(lambda, lambda_size, 1);
  const double c2 = sum(lambda, lambda_size, 2);
  const double c3 = sum(lambda, lambda_size, 3);
  const double c4 = sum(lambda, lambda_size, 4);
  double s1 = c3 / c2 / sqrt(c2);
  double s2 = c4 / c2 / c2;
  const double muQ = c1;
  const double sigmaQ = sqrt(2.0 * c2);
  const double tstar = (Q - muQ) / sigmaQ;

  double a;
  double delta;
  double l;
  if (s1 * s1 > s2) {
    a = 1 / (s1 - sqrt(s1 * s1 - s2));
    delta = (s1 * a - 1) * a * a;
    l = a * a - 2.0 * delta;
  } else {
    a = 1.0 / s1;
    delta = 0.0;
    l = c2 * c2 * c2 / c3 / c3;
  }
  const double muX = l + delta;
  const double sigmaX = sqrt(2) * a;

  int which = 1;
  double p;
  double q;
  double x = tstar * sigmaX + muX;
  double ncp = delta;
  int status;
  double bound;

  cdfchn(&which, &p, &q, &x, &l, &ncp, &status, &bound);
  if (status != 0) {
    return 1;
  }
  return q;
}

void MixtureChiSquare::dumpLambda() const {
  for (int i = 0; i < lambda_size; ++i) {
    fprintf(stderr, "lambda[%d] = %g\n", i, lambda[i]);
  }
}
