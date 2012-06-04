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

    double GetPvalue() const {return this->pvalue;};

  private:
    void splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol);
    double pvalue;
    LinearRegression lr;
};



#endif /* _LINEARREGRESSIONSCORETEST_H_ */
