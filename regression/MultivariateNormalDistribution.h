#ifndef _MULTIVARIATENORMAL_H_
#define _MULTIVARIATENORMAL_H_

#include <vector>
#include "MathMatrix.h"

//////////////////////////////////////////////////
// a C++ class for multivariate normal distribution
class MultivariateNormal{
 public:
  MultivariateNormal() {
    nu_ = 0;            // should be 0 for multivarate normal (degree of freedom)
    maxpts_ = 25000;    // default in mvtnorm: 25000
    abseps_ = 1e-6;     // default in mvtnorm: 0.001, we make it more stringent
    releps_ = 0;        // default in mvtnorm: 0
  }

  /**
   * @param n, number of variable
   * @param lower, lower limit
   * @param upper, upper limit
   * @param cor, correlation coefficient in row I column J stored
   *             in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
   * @param result, store results
   * @return 0 if succeed
   */
  int getBandProbFromCor(int n,
                         double* lower,
                         double* upper,
                         double* cov,
                         double* result);
  /**
   * @param n, number of variable
   * @param lower, lower limit
   * @param upper, upper limit
   * @param cov, covariance
   * @param result, store results
   * @return 0 if succeed
   */
  int getUpperFromCov(int n,
                      double* lower,
                      Matrix& cov,
                      double* result);

  int getUpperFromCov(double lower,
                      Matrix& cov,
                      double* result);


 private:
  int nu_ ;
  int maxpts_ ;         // default in mvtnorm: 25000
  double abseps_ ;      // default in mvtnorm: 0.001, we make it more stringent
  double releps_ ;      // default in mvtnorm: 0
  double error_;        // estimated abs. error with 99% confidence interval

  std::vector<double> lower;
  std::vector<double> upper;
  std::vector<int> infin;
  std::vector<double> delta;
  int inform;
};

#endif
