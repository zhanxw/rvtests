#ifndef _LINEARREGRESSIONSCORETEST_H_
#define _LINEARREGRESSIONSCORETEST_H_

#include <cmath>
#include "libsrc/MathMatrix.h"
#include "regression/LinearRegression.h"

class LinearRegressionScoreTest {
 public:
  LinearRegressionScoreTest();
  /**
   * @param colToTest: 0-based index, for that column of X will be tested
   */
  bool FitLinearModel(Matrix& X, Vector& y, int colToTest);

  bool FitNullModel(Matrix& Xnull, Vector& y);
  bool TestCovariate(Matrix& Xnull, Vector& y, Vector& Xcol);
  /**
   * Test H0: \beta = 0  (\beta is multiple dimension).
   * y ~ \beta * Xcol + \gamma * Xnull
   */
  bool TestCovariate(Matrix& Xnull, Vector& y, Matrix& Xcol);

  // fit y ~ 1 + beta*x  (no covariate)
  bool TestCovariate(Vector& x, Vector& y);
  /**
   * Test y ~ 1 + \beta * X (no covariate)
   */
  bool TestCovariate(const Matrix& x, const Vector& y);

  double GetStat() const { return this->stat; };
  double GetPvalue() const { return this->pvalue; };

  // U and V matrix are square matrice with dimension equalling to Xcol
  // e.g. when testing a single beta coefficient, U and V are 1 by 1 matrice.
  // e.g. when there is no covariates or intercept
  // U = x' * (y - y_hat)
  // V = x' * x * sigma2
  const Matrix& GetU() const { return this->Umatrix; };
  const Matrix& GetV() const { return this->Vmatrix; };
  const Matrix& GetBeta() const { return this->betaMatrix; }
  double GetSigma2() const { return this->lr.GetSigma2(); };
  const double GetSEBeta(int idx) const;

  // get estimates from null models
  Vector& GetNullCovEst() { return this->lr.GetCovEst(); };  // (X'X)^{-1} X'Y
  Matrix& GetNullCovB() { return this->lr.GetCovB(); };

 private:
  void splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol);
  double pvalue;
  double stat;
  Matrix Umatrix;
  Matrix Vmatrix;
  Matrix betaMatrix;
  LinearRegression lr;
};

#endif /* _LINEARREGRESSIONSCORETEST_H_ */
