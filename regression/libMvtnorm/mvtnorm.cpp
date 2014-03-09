#include "mvtnorm.h"

#include <stdlib.h>
#include <stdio.h>

// error message
const static char errorMessage0[] = "Normal Completion";                           // inform = 0
const static char errorMessage1[] = "Completion with error > abseps";              // inform = 1
const static char errorMessage2[] = "N greater 1000 or N < 1";                     // inform = 2
const static char errorMessage3[] ="Covariance matrix not positive semidefinite";  // inform = 3
const static char* errorMessage[4] = {errorMessage0, errorMessage1, errorMessage2, errorMessage3};


// infinity bounds
const static int INFIN_BOUND_NORMAL_NORMAL = 2;        // (..., ...)
const static int INFIN_BOUND_NORMAL_INFTY = 1;         // (..., inf)
const static int INFIN_BOUND_LOWER_NORMAL = 0;         // (-inf, ..)
const static int INFIN_BOUND_INFTY_INFTY = -1;  // (-inf, inf)

extern "C"{
  extern void mvtdst_(int* n,
                      int* nu,
                      double* lower,
                      double* upper,
                      int* infin,
                      double* correl,
                      double* delta,
                      int* maxpts,
                      double* abseps,
                      double* releps,
                      double* error,
                      double* value,
                      int* inform);
}
/**
 * @return <0 if anything goes wrong.
 */
int pmvnorm(int* n,
            int* nu,       // the number of degrees of freedom, if nu < 0, then use multivariate normal.
            double* lower,
            double* upper,
            int* infin,
            double* correl,
            double* delta, // non-central parameter
            int* maxpts,    // param
            double* abseps, // param
            double* releps, // param
            double* error,  // estimated abs. error. with 99% confidence interval
            double* value,     // results store here.
            int* inform)    // inform message goes here
{
  mvtdst_ (n, nu,
           lower, upper, infin, correl, delta,
           maxpts, abseps, releps,
           error, value, inform);
  // fprintf(stderr, "error = %g, value = %g, inform = %d\n", *error, *value, *inform);

  switch (*inform) {
    case 0:
      return *value;
    case 1:
    case 2:
    case 3:
      return -1.0;
  };

  return 0;
};

/**
 * @return CDF of multivariate normal P ( X < bound ) where X ~ MVN(0, correlationMatrix)
 */
double pmvnorm_P(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error)
{
  int nu_ = 0;
  int maxpts_ = 25000;     // default in mvtnorm: 25000
  double abseps_ = 1e-6;   // default in mvtnorm: 0.001, we make it more stringent
  double releps_ = 0;      // default in mvtnorm: 0

  double* lower = new double[n];
  int* infin = new int[n];
  double* delta = new double[n];

  int i = 0;
  for (i = 0; i < n; ++i) {
    infin[i] = 0; // (-inf, bound]
    lower[i] = 0.0;
    delta[i] = 0.0;
  }

  // return values
  double value_ = 0.0;
  int inform_ = 0.0;

   int ret = pmvnorm(&n, &nu_, lower, bound, infin, correlationMatrix, delta, &maxpts_, &abseps_, &releps_, error, &value_, &inform_);
  delete[] (lower);
  delete[] (infin);
  delete[] (delta);

  if (ret) {return -1.0;}
  return value_;
}

/**
 * @return (1 - CDF) of multivariate normal P ( X > bound ) where X ~ MVN(0, correlationMatrix)
 */
double pmvnorm_Q(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error)
{
  int nu_ = 0;
  int maxpts_ = 25000;
  double abseps_ = 1e-6;
  double releps_ = 0;


  double* upper = new double[n];
  int* infin = new int[n];
  double* delta = new double[n];

  int i = 0;
  for (i = 0; i < n; ++i) {
    infin[i] = 1; // (-inf, bound]
    upper[i] = 0.0;
    delta[i] = 0.0;
  }

  // return values
  double value_ = 0.0;
  int inform_ = 0.0;

  int ret = pmvnorm(&n, &nu_, bound, upper, infin, correlationMatrix, delta, &maxpts_, &abseps_, &releps_, error, &value_, &inform_);
  delete[] (upper);
  delete[] (infin);
  delete[] (delta);

  if (ret) {return -1.0;}
  return value_;
}

//////////////////////////////////////////////////
// a C++ class for multivariate normal distribution
int MultivariateNormal::getBandProbFromCov(int n,
                                           double* lower,
                                           double* upper,
                                           double* cov,
                                           double* result) {
  infin.resize(n);
  std::fill(infin.begin(), infin.end(), INFIN_BOUND_NORMAL_NORMAL);
  delta.resize(n);
  std::fill(delta.begin(), delta.end(), 0);

  return pmvnorm(&n,
                 &nu_,
                 lower,
                 upper,
                 infin.data(),
                 cov,
                 delta.data(),
                 &maxpts_,
                 &abseps_,
                 &releps_,
                 &error_,
                 result,
                 &inform);
};
