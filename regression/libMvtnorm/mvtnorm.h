#ifndef _MVT_H_
#define _MVTNORM_H_

#include <vector>

/**
 * @return 0 if succeed
 */
int pmvnorm(int* n,
            int* nu,
            double* lower,
            double* upper,
            int* infin,      //   if INFIN(I) < 0, Ith limits are (-infinity, infinity);
            //   if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
            //   if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
            //   if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
            double* correl,  //   the correlation coefficient in row I column J of the
            //     correlation matrixshould be stored in
            //    CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
            double* delta,   // non-central parameter
            int* maxpts,     // param
            double* abseps,  // param
            double* releps,  // param
            double* error,   // estimated abs. error. with 99% confidence interval
            double* value,      // results store here.
            int* inform);    // inform message goes here

/**
 * Quick way to get pvalue
 * @return p-value, or -1.0 if error happens
 */
double pmvnorm_P(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error);

/**
 * Quick way to get pvalue
 * @return p-value, or -1.0 if error happens
 */
double pmvnorm_Q(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error);


class MultivariateNormal{
 public:
  MultivariateNormal() {
    nu_ = 0;            //
    maxpts_ = 25000;    // default in mvtnorm: 25000
    abseps_ = 1e-6;     // default in mvtnorm: 0.001, we make it more stringent
    releps_ = 0;        // default in mvtnorm: 0
  }
  int getBandProbFromCov(int n,
                         double* lower,
                         double* upper,
                         double* cov,
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

#endif /* _MVTNORM_H_ */
