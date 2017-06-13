#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>

#include "GSLMinimizer.h"

Minimizer::Minimizer()
    : finalX(NAN), finalY(NAN), epsabs(0.001), epsrel(0.0), maxIter(100) {
  T = gsl_min_fminimizer_brent;
  s = gsl_min_fminimizer_alloc(T);
}

Minimizer::~Minimizer() { gsl_min_fminimizer_free(s); }

int Minimizer::minimize(gsl_function F, double startValue, double lowerBound,
                        double upperBound) {
  gsl_error_handler_t* oldHandler = gsl_set_error_handler_off();
  int iter = 0;
  double a, b;

  status = gsl_min_fminimizer_set(s, &F, startValue, lowerBound, upperBound);
  if (status != GSL_SUCCESS) {
#ifndef NDEBUG
    fprintf(stderr, "Minimizer failed due to: %s\n", gsl_strerror(status));
#endif
    goto errorLabel;
  }

  do {
    iter++;
    status = gsl_min_fminimizer_iterate(s);
    if (status == GSL_EBADFUNC || status == GSL_FAILURE) {
      goto errorLabel;
    }

    finalX = gsl_min_fminimizer_x_minimum(s);
    a = gsl_min_fminimizer_x_lower(s);
    b = gsl_min_fminimizer_x_upper(s);

    // stopping rule
    // see
    // http://www.gnu.org/software/gsl/manual/html_node/Minimization-Stopping-Parameters.html
    status = gsl_min_test_interval(a, b, epsabs, epsrel);

    if (status == GSL_SUCCESS) {
      // printf ("Converged:\n");
      finalY = gsl_min_fminimizer_f_minimum(s);
      goto successLabel;
    }
    // printf("%5d\ta=%g\tb=%g\tfinalX=%g\n", iter, a, b, finalX);
    // printf ("%5d .7f, [.7f] "
    //         "%.7f %+.7f %.7f\n",
    //         iter, a, b,
    //         m, m - m_expected, b - a);
  } while (status == GSL_CONTINUE && iter < this->maxIter);

successLabel:
  gsl_set_error_handler(oldHandler);
  return 0;
errorLabel:
  gsl_set_error_handler(oldHandler);
  return -1;
}
