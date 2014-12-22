#include "mvtnorm.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <vector>
 const int MvtNorm::INFIN_BOUND_NORMAL_NORMAL = 2;        // (..., ...)
 const int MvtNorm::INFIN_BOUND_NORMAL_INFTY = 1;         // (..., inf)
 const int MvtNorm::INFIN_BOUND_LOWER_NORMAL = 0;         // (-inf, ..)
 const int MvtNorm::INFIN_BOUND_INFTY_INFTY = -1;         // (-inf, inf)

const char* MvtNorm::informMessage[] = {
  "Normal Completion",                           // inform = 0
  "Completion with error > abseps",              // inform = 1
  "N greater 1000 or N < 1",                     // inform = 2
  "Covariance matrix not positive semidefinite"  // inform = 3
};

// // error message
// const static char errorMessage0[] = "Normal Completion";                           // inform = 0
// const static char errorMessage1[] = "Completion with error > abseps";              // inform = 1
// const static char errorMessage2[] = "N greater 1000 or N < 1";                     // inform = 2
// const static char errorMessage3[] ="Covariance matrix not positive semidefinite";  // inform = 3
// const static char* errorMessage[4] = {errorMessage0, errorMessage1, errorMessage2, errorMessage3};

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

MvtNorm::MvtNorm():
    nu(0),
    maxpts(25000),     // default in mvtnorm: 25000
    abseps(0.001),   // default in mvtnorm: 0.001
    releps(0.)      // default in mvtnorm: 0
{
}

MvtNorm::MvtNorm(int nu, int maxpts, double abseps, double releps):
    nu(0),
    maxpts(25000),     // default in mvtnorm: 25000
    abseps(0.001),   // default in mvtnorm: 0.001
    releps(0.)      // default in mvtnorm: 0
{
}

/**
 * Quick way to get pvalue
 * @return CDF of multivariate normal P ( lower < X < upper ) where X ~ MVN(0, correlationMatrix)
 * @return p-value, or -1.0 if error happens
 */
double MvtNorm::compute_Band(int n,
                             double* lower,
                             double* upper,
                             double* correlationMatrix) // (2,1), (3,1), (3,2) .....
{
  infin.resize(n);
  std::fill(infin.begin(), infin.end(), MvtNorm::INFIN_BOUND_NORMAL_NORMAL);

  int ret = compute(&n, lower, upper, infin.data(), correlationMatrix);

  if (ret) {return -1.0;}
  return prob;
}


/**
 * @return CDF of multivariate normal P ( X < bound ) where X ~ MVN(0, correlationMatrix)
 */
double MvtNorm::compute_P(int n,
                          double* bound,
                          double* correlationMatrix) // (2,1), (3,1), (3,2) .....
{
  lower.resize(n);
  infin.resize(n);
  std::fill(infin.begin(), infin.end(), MvtNorm::INFIN_BOUND_LOWER_NORMAL);

  int ret = compute(&n, lower.data(), bound, infin.data(), correlationMatrix);

  if (ret) {return -1.0;}
  return prob;
}

/**
 * @return (1 - CDF) of multivariate normal P ( X > bound ) where X ~ MVN(0, correlationMatrix)
 */
double MvtNorm::compute_Q(int n,
                          double* bound,
                          double* correlationMatrix) // (2,1), (3,1), (3,2) .....
{
  upper.resize(n);
  infin.resize(n);
  std::fill(infin.begin(), infin.end(), MvtNorm::INFIN_BOUND_NORMAL_INFTY);

  int ret = compute(&n, bound, upper.data(), infin.data(), correlationMatrix);

  if (ret) {return -1.0;}
  return prob;
}


/**
 * @return 0 if succeed
 */
int MvtNorm::compute(int* n,
                     double* lower,
                     double* upper,
                     int* infin,
                     double* correl)
{
#ifdef DEBUG
  const int N = *n;
  fprintf(stderr, "n = %d, ", N);
  fprintf(stderr, "lower = ");
  for (int i = 0; i < N; ++i) {
    fprintf(stderr, "%g, ", lower[i]);
  }
  fprintf(stderr, "upper = ");
  for (int i = 0; i < N; ++i) {
    fprintf(stderr, "%g, ", upper[i]);
  }
  fprintf(stderr, "correl = ");
  for (int i = 0; i < N*(N-1)/2; ++i) {
    fprintf(stderr, "%g, ", correl[i]);
  }
#endif

  delta.resize(*n);
  std::fill(delta.begin(), delta.end(), 0.);
  
  mvtdst_ (n, &this->nu,
           lower, upper, infin, correl, delta.data(),
           &this->maxpts, &this->abseps, &this->releps,
           &this->error, &this->prob, &this->inform);

#ifdef DEBUG
  fprintf(stderr, "error = %g, value = %g, inform = %d\n", error, prob, inform);
#endif

  return (inform) ;
}

const char* MvtNorm::getError(int i) {
  const int n = sizeof(MvtNorm::informMessage)/sizeof(MvtNorm::informMessage[0]);
  if (i <= 0 || i >= n) return NULL;
  return MvtNorm::informMessage[i];
}
