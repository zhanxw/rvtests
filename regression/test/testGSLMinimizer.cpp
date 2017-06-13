#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "regression/GSLMinimizer.h"
#include "third/gsl/include/gsl/gsl_errno.h"

static double func1(double x, void* param) { return x * x; }

static double func2(double x, void* param) {
  // this example is taken from:
  // https://stackoverflow.com/questions/18807512/is-this-a-bug-in-gsls-minimum-finding-routine

  // direct optimization will make it fail
  double beta = x;
  return -(
      ((6160558822864 * (exp(4 * beta)) + 523830424923 * (exp(3 * beta)) +
        1415357447750 * (exp(5 * beta)) + 7106224104 * (exp(6 * beta))) /
       (385034926429 * (exp(4 * beta)) + 58203380547 * (exp(3 * beta)) +
        56614297910 * (exp(5 * beta)) + 197395114 * (exp(6 * beta)))) -
      ((1540139705716 * (exp(4 * beta)) + 174610141641 * (exp(3 * beta)) +
        283071489550 * (exp(5 * beta)) + 1184370684 * (exp(6 * beta))) *
       (1540139705716 * (exp(4 * beta)) + 174610141641 * (exp(3 * beta)) +
        283071489550 * (exp(5 * beta)) + 1184370684 * (exp(6 * beta))) /
       pow((385034926429 * (exp(4 * beta)) + 58203380547 * (exp(3 * beta)) +
            56614297910 * (exp(5 * beta)) + 197395114 * (exp(6 * beta))),
           2)));
}

int main() {
  {
    Minimizer m;
    gsl_function F;
    F.function = func1;
    F.params = NULL;
    double start = 0.5;
    double lb = -5;
    double ub = 5;
    if (m.minimize(F, start, lb, ub)) {
      fprintf(stderr, "Minimizer failed");
      assert(false);
    } else {
      fprintf(stderr, "Minimizer succeed, x = %g\n", m.getX());
    }
  }
  {
    Minimizer m;
    gsl_function F;
    F.function = func2;
    F.params = NULL;
    double start = 0.0;
    double lb = -6;
    double ub = 6;

    // gsl_set_error_handler_off();

    if (m.minimize(F, start, lb, ub)) {
      fprintf(stderr, "Minimizer failed\n");
    } else {
      fprintf(stderr, "Minimizer succeed, x = %g\n", m.getX());
      assert(false);
    }
  }

  return 0;
}
