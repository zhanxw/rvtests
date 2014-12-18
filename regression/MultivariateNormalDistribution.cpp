#include "MultivariateNormalDistribution.h"

#include <math.h>

#include "libMvtnorm/mvtnorm.h"
#include "MathMatrix.h"

int MultivariateNormal::getBandProbFromCor(int n,
                                           double* lower,
                                           double* upper,
                                           double* cor,
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
                 cor,
                 delta.data(),
                 &maxpts_,
                 &abseps_,
                 &releps_,
                 &error_,
                 result,
                 &inform);
};

int MultivariateNormal::getUpperFromCov(int n,
                                        double* lower,
                                        Matrix& cov,
                                        double* result) {
  std::vector<double> v(n, 0.);
  for (int i = 0; i < n; ++i) {
    v[i] = sqrt((cov)[i][i]);
  }
  
  std::vector<double> cor(n*(n-1)/2, 0.);
  int k = 0;
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      cor[k++] = (cov)[i][j] / v[i] / v[j];
    }
  }
  assert(k == (n * (n - 1) / 2));

  return pmvnorm_Q(n, lower, cor.data(), &error_)  ;
}


int MultivariateNormal::getUpperFromCov(double lower,
                                        Matrix& cov,
                                        double* result) {
  if (cov.rows != cov.cols) {
    return -1;
  }
  
  Matrix& m = cov;
  int n = m.rows;

  std::vector<double> lowerBound(n, lower);
  return this->getUpperFromCov(n,
                               lowerBound.data(),
                               cov,
                               result);
}

