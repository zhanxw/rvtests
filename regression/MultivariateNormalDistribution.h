#ifndef _MULTIVARIATENORMALDISTRIBUTION_H_
#define _MULTIVARIATENORMALDISTRIBUTION_H_

#include <vector>
#include "MathMatrix.h"
#include "libMvtnorm/mvtnorm.h"

//////////////////////////////////////////////////
// a C++ class for multivariate normal distribution
class MultivariateNormalDistribution{
 public:
  double getAbsEps() const {
    return mvn.getAbsEps();
  }
  //////////////////////////////////////////////////
  // Get band probability
  //////////////////////////////////////////////////  
  /**
   * @param n, number of variable
   * @param lower, lower limit
   * @param upper, upper limit
   * @param cor, correlation coefficient in row I column J stored
   *             in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
   * @param result, store results
   * @return 0 if succeed, 1 if succeed with bigger error, otherwise means fail
   */
  int getBandProbFromCor(int n,
                         double* lower,
                         double* upper,
                         double* cor,
                         double* result);
  int getBandProbFromCor(int n,
                         double* lower,
                         double* upper,
                         Matrix& cor,
                         double* result);
  int getBandProbFromCov(int n,
                         double* lower,
                         double* upper,
                         Vector& mean,
                         Matrix& cov,
                         double* result);

  //////////////////////////////////////////////////
  // Get upper probability
  //////////////////////////////////////////////////    
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
                      Vector& mean,
                      Matrix& cov,
                      double* result);
  // similar to above, except assuming the means are zeros  
  int getUpperFromCov(int n,
                      double* lower,
                      Matrix& cov,
                      double* result);
  // similar to above, except assuming the means are zeros  
  int getUpperFromCov(double lower,
                      Matrix& cov,
                      double* result);

  //////////////////////////////////////////////////
  // Get lower probability
  //////////////////////////////////////////////////  
  int getLowerFromCov(int n,
                      double* upper,
                      Vector& mean,
                      Matrix& cov,
                      double* result);
  int getLowerFromCov(int n,
                      double* upper,
                      Matrix& cov,
                      double* result);
  int getLowerFromCov(double upper,
                      Matrix& cov,
                      double* result);

 private:
  /**
   * Convert a covariance matrix @param m to correlation matrix,
   * and store it to a 1-d vector @param out.
   */
  void toCor(Matrix& m, std::vector<double>* out);

 private:
  MvtNorm mvn;
};

#endif
