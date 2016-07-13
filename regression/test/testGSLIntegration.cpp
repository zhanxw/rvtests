#include <stdio.h>
#include "GSLIntegration.h"

double f(double x, void* param) {
  if (x > 0 ) {
    return (exp(-x));
  }
  return (exp(x));
}

int main(int argc, char *argv[])
{
  gsl_function F;
  F.function = &f;
  
  Integration integral;
  if (integral.integrate(F)) {
    return 1;
  }

  fprintf(stderr, "expected = %g\n", 2.0);
  fprintf(stderr, "actual = %g\n", integral.getResult());
  
  return 0;
}

