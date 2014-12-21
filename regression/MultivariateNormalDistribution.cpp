#include "MultivariateNormalDistribution.h"

#include <math.h>
#include <gsl/gsl_cdf.h>


#include "MathMatrix.h"

//////////////////////////////////////////////////
// Get band probability
//////////////////////////////////////////////////
int MultivariateNormalDistribution::getBandProbFromCor(int n,
                                                       double* lower,
                                                       double* upper,
                                                       double* cor,
                                                       double* result) {
  if (n == 1) {
    *result =  (gsl_cdf_ugaussian_P(*upper) -
                gsl_cdf_ugaussian_P(*lower));
    return 0;
  }

  *result = mvn.compute_Band(n, lower, upper, cor);
  return (*result >= 0. ? 0 : -1);
}

int MultivariateNormalDistribution::getBandProbFromCor(int n,
                                                       double* lower,
                                                       double* upper,
                                                       Matrix& cor,
                                                       double* result) {
  // skip checking diagnol elements equaling to one.

  std::vector<double> tmp;
  toCor(cor, &tmp);
  return getBandProbFromCor(n, lower, upper, tmp.data(), result);
}

int MultivariateNormalDistribution::getBandProbFromCov(int n,
                                                       double* lower,
                                                       double* upper,
                                                       Vector& mean,
                                                       Matrix& cov,
                                                       double* result) {
  // FIX here and other place, divide diag(v)
  for (int i = 0; i < n; ++i) {
    (lower)[i] -= mean[i]/sqrt(cov[i][i]);
    (upper)[i] -= mean[i]/sqrt(cov[i][i]);
  }
  std::vector<double> cor;
  toCor(cov, &cor);
  return getBandProbFromCor(n, lower, upper, cor.data(), result);
}

//////////////////////////////////////////////////
// Get upper probability
//////////////////////////////////////////////////
int MultivariateNormalDistribution::getUpperFromCov(int n,
                                                    double* lower,
                                                    Matrix& cov,
                                                    double* result) {
  if (n == 1) {
    *result =  (gsl_cdf_gaussian_Q(*lower, sqrt(cov[0][0])));
    return 0;
  }

  for (int i = 0; i < n; ++i) {
    lower[i] /= sqrt(cov[i][i]);
  }
  std::vector<double> cor;
  toCor(cov, &cor);

  *result = mvn.compute_Q(n, lower, cor.data());
  return (*result >= 0. ? 0 : -1);
}

int MultivariateNormalDistribution::getUpperFromCov(int n,
                                                    double* lower,
                                                    Vector& mean,
                                                    Matrix& cov,
                                                    double* result) {
  for (int i = 0; i < n; ++i) {
    lower[i] -= mean[i];
  }
  return getUpperFromCov(n, lower, cov, result);
}

int MultivariateNormalDistribution::getUpperFromCov(double lower,
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
//////////////////////////////////////////////////
// Get lower probability
//////////////////////////////////////////////////
int MultivariateNormalDistribution::getLowerFromCov(int n,
                                                    double* upper,
                                                    Matrix& cov,
                                                    double* result) {
  if (n == 1) {
    *result = (gsl_cdf_gaussian_P(*upper, sqrt(cov[0][0])));
    return 0;
  }

  for (int i = 0; i < n; ++i) {
    upper[i] /= sqrt(cov[i][i]);
  }
  std::vector<double> cor;
  toCor(cov, &cor);

  *result = mvn.compute_P(n, upper, cor.data());
  return (*result >= 0. ? 0 : -1);
}

int MultivariateNormalDistribution::getLowerFromCov(int n,
                                                    double* upper,
                                                    Vector& mean,
                                                    Matrix& cov,
                                                    double* result) {
  for (int i = 0; i < n; ++i) {
    upper[i] -= mean[i];
  }
  return getLowerFromCov(n, upper, cov, result);
}

int MultivariateNormalDistribution::getLowerFromCov(double upper,
                                                    Matrix& cov,
                                                    double* result) {
  if (cov.rows != cov.cols) {
    return -1;
  }

  Matrix& m = cov;
  int n = m.rows;

  std::vector<double> upperBound(n, upper);
  return this->getLowerFromCov(n,
                               upperBound.data(),
                               cov,
                               result);
}

//////////////////////////////////////////////////
// Utility functions
//////////////////////////////////////////////////
void MultivariateNormalDistribution::toCor(Matrix& m, std::vector<double>* out) {
  if (m.rows != m.cols) return;
  if (m.rows == 0) return;

  const int n = m.rows;

  std::vector<double> v(n, 0.);
  for (int i = 0; i < n; ++i) {
    v[i] = sqrt((m)[i][i]);
  }

  std::vector<double>& cor = (*out);
  cor.resize(n*(n-1)/2, 0.);
  int k = 0;
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      cor[k++] = (m)[i][j] / v[i] / v[j];
      // fprintf(stderr, "%g\n", (m)[i][j] / v[i] / v[j]);
    }
  }
  assert(k == (n * (n - 1) / 2));
}
