#ifndef GSLINTEGRATION_H
#define GSLINTEGRATION_H

#include <gsl/gsl_integration.h>

class Integration
{
 public:
  Integration();
  virtual ~Integration();

  // integrate from -\inf to \inf
  int integrate(gsl_function F);
  double getResult() const {
    return this->result;
  }
  // get an estimate of the absolute error
  double getAbsError() const {
    return this->abserr;
  }
 private:
  // GSL stuffs
  double epsabs;
  double epsrel;
  size_t limit;
  gsl_integration_workspace * workspace;
  double result;
  double abserr;
};

#endif /* GSLINTEGRATION_H */
