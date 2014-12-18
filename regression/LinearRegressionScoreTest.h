#ifndef _LINEARREGRESSIONSCORETEST_H_
#define _LINEARREGRESSIONSCORETEST_H_

#include "MathMatrix.h"
#include "MathCholesky.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

#include "LinearRegression.h"

class LinearRegressionScoreTest{
 public:
  LinearRegressionScoreTest();
  /**
   * @param colToTest: 0-based index, for that column of X will be tested
   */
  bool FitLinearModel(Matrix &X, Vector &y, int colToTest);

  bool FitNullModel(Matrix& Xnull, Vector& y);
  bool TestCovariate(Matrix& Xnull, Vector& y, Vector& Xcol);
  /**
   * Test H0: \beta = 0  (\beta is multiple dimension).
   * y ~ \beta * Xcol + \gamma * Xnull
   */
  bool TestCovariate(Matrix& Xnull, Vector& y, Matrix& Xcol);

  // fit y~1+ beta*x  (no covariate)
  bool TestCovariate(Vector& x, Vector& y);
  /**
   * Test y~ 1 + \beta * X (no covariate)
   */
  bool TestCovariate(Matrix& x, Vector& y);

  double GetStat() const {return this->stat;};
  double GetPvalue() const {return this->pvalue;};

  // U and V matrix are square matrice with dimension equalling to Xcol
  // e.g. when testing a single beta coefficient, U and V are 1 by 1 matrice.
  const Matrix& GetU() const {return this->Umatrix;};
  const Matrix& GetV() const {return this->Vmatrix;};
  const Matrix& GetBeta() const { return this->betaMatrix;}
  const double GetSigma2() const { return this->lr.GetSigma2();};
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
