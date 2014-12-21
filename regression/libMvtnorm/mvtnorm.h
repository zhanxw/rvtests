#ifndef _MVTNORM_H_
#define _MVTNORM_H_

#include <vector>

/**
 * a low-level class to wrap mvtnorm distributions
 */
class MvtNorm{
 public:
  MvtNorm();
  MvtNorm(int nu, int maxpts, double abseps, double releps);
  /**
   * Quick way to get prob
   * @return CDF of multivariate normal P ( lower < X < upper ) where X ~ MVN(0, correlationMatrix)
   * @return p-value, or -1.0 if error happens
   */
  double compute_Band(int n,
                      double* lower,
                      double* upper,
                      double* correlationMatrix); // (2,1), (3,1), (3,2) .....

  /**
   * Quick way to get prob
   * @return CDF of multivariate normal P ( X < bound ) where X ~ MVN(0, correlationMatrix)
   * @return p-value, or -1.0 if error happens
   */
  double compute_P(int n,
                   double* bound,
                   double* correlationMatrix // (2,1), (3,1), (3,2) .....
                   );

  /**
   * Quick way to get prob
   * @return (1 - CDF) of multivariate normal P ( X > bound ) where X ~ MVN(0, correlationMatrix)
   * @return p-value, or -1.0 if error happens
   */
  double compute_Q(int n,
                   double* bound,
                   double* correlationMatrix // (2,1), (3,1), (3,2) .....
                   );

  const char* getError(int i);
  double getAbsEps() const {
    return this->abseps;
  }
 private:
  /**
   * @return <0 if anything goes wrong.
   */
  int compute(int* n,
              double* lower,
              double* upper,
              int* infin,
              double* correl);
 public:
  // infinity bounds
  static const int INFIN_BOUND_NORMAL_NORMAL;
  static const int INFIN_BOUND_NORMAL_INFTY ;
  static const int INFIN_BOUND_LOWER_NORMAL; 
  static const int INFIN_BOUND_INFTY_INFTY;
 private:
  // error message
  // static const char errorMessage0[] = "Normal Completion";                           // inform = 0
  // static const char errorMessage1[] = "Completion with error > abseps";              // inform = 1
  // static const char errorMessage2[] = "N greater 1000 or N < 1";                     // inform = 2
  // static const char errorMessage3[] ="Covariance matrix not positive semidefinite";  // inform = 3
  static const char* informMessage[4];

  // calculation results
  double prob;
  double error;            // error of the calculated prob
  int inform;           // calculation status

  // temporary values
  std::vector<double> upper;
  std::vector<double> lower;
  std::vector<int> infin;
  std::vector<double> delta;

  // parameters
  int nu;         // the number of degrees of freedom, if nu < 0, then use multivariate normal.
  int maxpts;     // default in mvtnorm: 25000
  double abseps;   // default in mvtnorm: 0.001
  double releps;      // default in mvtnorm: 0

};

#endif /* _MVTNORM_H_ */
