#include <math.h>
#include <stdio.h>

#include <gsl/gsl_cdf.h>

#include "IO.h"

#include "MathMatrix.h"
#include "MathVector.h"
#include "MatrixIO.h"
#include "MatrixOperation.h"

#include "MultivariateNormalDistribution.h"

bool equal(double a, double b, double eps) {
  return (fabs(a - b) <= eps);
}

int main(int argc, char *argv[])
{
  MultivariateNormalDistribution mvn;
  const double eps = mvn.getAbsEps();
  
  // repeat process by pmvnorm() in libmvtnorm R package
  {
    /*
       n <- 5
       mean <- rep(0, 5)
       lower <- rep(-1, 5)
       upper <- rep(3, 5)
       corr <- diag(5)
       corr[lower.tri(corr)] <- 0.5
       corr[upper.tri(corr)] <- 0.5
       prob <- pmvnorm(lower, upper, mean, corr)
       print(prob)
    */
    
    Matrix cor;
    cor.Dimension(5, 5);
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < 5; ++j) {
        if (i == j) {
          cor[i][j] = 1.;
        } else {
          cor[i][j] = .5;
        }
      }
    }
    std::vector<double> lower(5, -1.);
    std::vector<double> upper(5, 3.);

    double res = -1.;
    int ret = mvn.getBandProbFromCor(5,
                                     lower.data(),
                                     upper.data(),
                                     cor,
                                     &res);
    assert(ret == 0);
    assert(equal(res, 0.5800051, eps));
  }

  {
    
    /*

      stopifnot(pmvnorm(lower=-Inf, upper=3, mean=0, sigma=1) == pnorm(3))

    */
    Matrix v;
    v.Dimension(1, 1);
    v[0][0] = 1.;
    double res;
    int ret = mvn.getLowerFromCov(3., v, &res);
    assert(ret == 0);
    assert(equal(res, gsl_cdf_ugaussian_P(3.), eps) );
  }

  {
    /*
      a <- pmvnorm(lower=-Inf,upper=c(.3,.5),mean=c(2,4),diag(2))

      stopifnot(round(a,16) == round(prod(pnorm(c(.3,.5),c(2,4))),16))

     */
    std::vector<double> upper;
    upper.resize(2);
    upper[0] = .3;
    upper[1] = .5;
    Vector mean;
    mean.Dimension(2);
    mean[0] = 2.;
    mean[1] = 4.;
    Matrix v;
    v.Dimension(2, 2);
    v[0][0] = v[1][1] = 1.;
    v[1][0] = v[0][1] = 0.;

    double res;
    int ret = mvn.getLowerFromCov(2, upper.data(), mean, v, &res);
    assert(ret == 0);
    double truth = gsl_cdf_ugaussian_P(.3 - 2.) *
                    gsl_cdf_ugaussian_P(.5 - 4.);
    assert(equal(res, truth, eps));
    
  }

  {
    /*

      a <- pmvnorm(lower=-Inf,upper=c(.3,.5,1),mean=c(2,4,1),diag(3))

      stopifnot(round(a,16) == round(prod(pnorm(c(.3,.5,1),c(2,4,1))),16))
    */
    std::vector<double> upper;
    upper.resize(3);
    upper[0] = .3;
    upper[1] = .5;
    upper[2] = .1;
    Vector mean;
    mean.Dimension(3);
    mean[0] = 2.;
    mean[1] = 4.;
    mean[2] = 1.;
    Matrix v;
    v.Dimension(3, 3);
    v.Zero();
    v[0][0] = v[1][1] = v[2][2] = 1.;
    
    double res;
    int ret = mvn.getLowerFromCov(3, upper.data(), mean, v, &res);
    assert(ret == 0);
    double truth = gsl_cdf_ugaussian_P(.3 - 2.) *
                   gsl_cdf_ugaussian_P(.5 - 4.) *
                   gsl_cdf_ugaussian_P(1 - 1.);
    assert(ret == 0);
    assert(equal(res, truth, eps));
  }
  {
    /*
      # Example from R News paper (original by Genz, 1992):

      m <- 3
      sigma <- diag(3)
      sigma[2,1] <- 3/5
      sigma[3,1] <- 1/3
      sigma[3,2] <- 11/15
      pmvnorm(lower=rep(-Inf, m), upper=c(1,4,2), mean=rep(0, m), corr=sigma)
    */
    std::vector<double> upper;
    upper.resize(3);
    upper[0] = 1;
    upper[1] = 4;
    upper[2] = 2;
    Matrix v;
    v.Dimension(3, 3);
    v.Zero();
    v[1][0] = v[0][1] = 3./5;
    v[2][0] = v[0][2] = 1./3;
    v[2][1] = v[1][2] = 11./15;
    v[0][0] = v[1][1] = v[2][2] = 1.0;
    
    double res;
    int ret = mvn.getLowerFromCov(3, upper.data(), v, &res);
    assert(ret == 0);
    double truth = 0.8279851;
    assert(equal(res, truth, eps));
  }
  {
    /*
      # Correlation and Covariance

      a <- pmvnorm(lower=-Inf, upper=c(2,2), sigma = diag(2)*2)
      b <- pmvnorm(lower=-Inf, upper=c(2,2)/sqrt(2), corr=diag(2))
      stopifnot(all.equal(round(a,5) , round(b, 5)))
      
     */
    Matrix v;
    v.Dimension(2, 2);
    v.Zero();
    v[0][0] = v[1][1] = 2.;
    double upper = 2.;
    double res1;
    int ret = mvn.getLowerFromCov(2., v, &res1);
    assert(ret == 0);

    double res2;
    v[0][0] = v[1][1] = 1.0;
    ret = mvn.getLowerFromCov(2. / sqrt(2.), v, &res2);
    assert(ret == 0);
    assert(equal(res1, res2, eps) );
    
  }
  return 0;
};

