#include <math.h>
#include <stdio.h>

#include "IO.h"

#include "LogisticRegressionVT.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "MatrixIO.h"
#include "MatrixOperation.h"

int main(int argc, char *argv[]) {
  Vector Y;
  Matrix X;
  Matrix Cov;
  LoadVector("input.logistic.mvt.y", Y);
  LoadMatrix("input.logistic.mvt.x", X);
  LoadMatrix("input.logistic.mvt.cov", Cov);

  {
    Matrix x;
    Vector y;
    x = X;
    y = Y;

    LogisticRegressionVT logistic;
    if (!logistic.FitNullModel(Cov, Y)) {
      fprintf(stderr, "Fitting failed! - step 1!\n");
      return -1;
    }
    if (!logistic.TestCovariate(Cov, Y, X)) {
      fprintf(stderr, "Fitting failed - step 2!\n");
      return -1;
    }

    dumpToFile(logistic.GetU(), stdout);
    dumpToFile(logistic.GetV(), stdout);
    dumpToFile(logistic.GetT(), stdout);
    dumpToFile(logistic.GetCov(), stdout);
    fprintf(stdout, "%g\t0\n", logistic.GetPvalue(), 0);
  }
  return 0;
};
