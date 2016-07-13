#include <stdio.h>
#include <stdlib.h>
#include "MixtureChiSquare.h"

int main(int argc, char *argv[])
{
  {
    MixtureChiSquare mcs;
    mcs.addLambda(1.0);
    mcs.addLambda(2.0);
    mcs.addLambda(3.0);
    double p1 = mcs.getPvalue(4.0);
    double p2 = mcs.getLiuPvalue(4.0);
    const double p1Correct = 0.5510657;
    const double p2Correct = 0.5399171;
    fprintf(stderr, "%g\t%g\n", p1, p1Correct);
    fprintf(stderr, "%g\t%g\n", p2, p2Correct);
  }
  {
    MixtureChiSquare mcs;
    mcs.addLambda(1.0);
    mcs.addLambda(1.0);
    mcs.addLambda(1.0);
    double p1 = mcs.getPvalue(30.0);
    double p2 = mcs.getLiuPvalue(30.0);
    const double p1Correct = 2.557879e-06;
    const double p2Correct = 1.380057e-06;
    fprintf(stderr, "%g\t%g\n", p1, p1Correct);
    fprintf(stderr, "%g\t%g\n", p2, p2Correct);
  }
  {
    MixtureChiSquare mcs;
    mcs.addLambda(1.0);
    mcs.addLambda(1.0);
    mcs.addLambda(1.0);
    double p1 = mcs.getPvalue(50.0);
    double p2 = mcs.getLiuPvalue(50.0);
    const double p1Correct = 0.0;
    const double p2Correct = 7.989179e-11;
    fprintf(stderr, "%g\t%g\n", p1, p1Correct);
    fprintf(stderr, "%g\t%g\n", p2, p2Correct);
  }
  return 0;
}
