#include <math.h>

#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"
#include "FirthRegression.h"
#include "MatrixIO.h"

int main(int argc, char *argv[])
{
  Vector Y;
  Matrix X;

  LoadVector("input.firth.y", Y);
  LoadMatrix("input.firth.x", X);

  {
    Matrix x;
    Vector y;
    x = X;
    y = Y;

    FirthRegression firth;
    if ( firth.Fit(x, y) == false) {
      fprintf(stderr, "Fitting failed!\n");
      return -1;
    }
    Vector& beta  = firth.GetCovEst();
    Matrix& cov = firth.GetCovB();
    int n = beta.Length();
    for (int i = 0; i < n; ++i) {
      printf("%g\t%g\n", beta[i], sqrt(cov[i][i]));
    }
    // fprintf(stdout, "score_p\t");
    // double score_p =lrst.GetPvalue();
    // Print(score_p);
    // fputc('\n', stdout);
  }
  return 0;
};

