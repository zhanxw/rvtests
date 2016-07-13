#include <math.h>
#include <stdio.h>

#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"
#include "LinearRegressionVT.h"
#include "MatrixIO.h"
#include "MatrixOperation.h"

int main(int argc, char *argv[])
{
  Vector Y;
  Matrix X;
  Matrix Cov;
  LoadVector("input.linear.mvt.y", Y);
  LoadMatrix("input.linear.mvt.x", X);
  LoadMatrix("input.linear.mvt.cov", Cov);  

  {
    Matrix x;
    Vector y;
    x = X;
    y = Y;

    LinearRegressionVT linear;
    if (!linear.FitNullModel(Cov, Y)){
      fprintf(stderr, "Fitting failed! - step 1!\n");
      return -1;
    }
    if (!linear.TestCovariate(Cov, Y, X)){
      fprintf(stderr, "Fitting failed - step 2!\n");
      return -1;
    }

    dumpToFile(linear.GetU(), stdout);
    dumpToFile(linear.GetV(), stdout);
    dumpToFile(linear.GetT(), stdout);
    dumpToFile(linear.GetCov(), stdout);
    fprintf(stdout, "%g\t0\n", linear.GetPvalue(), 0);
  }
  return 0;
};

