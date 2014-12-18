#ifndef _MVT_H_
#define _MVTNORM_H_

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
 * @return CDF of multivariate normal P ( X < bound ) where X ~ MVN(0, correlationMatrix)
 * @return p-value, or -1.0 if error happens
 */
double pmvnorm_P(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error);

/**
 * Quick way to get pvalue
  * @return (1 - CDF) of multivariate normal P ( X > bound ) where X ~ MVN(0, correlationMatrix)
  * @return p-value, or -1.0 if error happens
 */
double pmvnorm_Q(int n,
                 double* bound,
                 double* correlationMatrix, // (2,1), (3,1), (3,2) .....
                 double* error);


#endif /* _MVTNORM_H_ */
