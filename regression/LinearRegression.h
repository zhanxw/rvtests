#ifndef _LINEARREGRESSION_H_
#define _LINEARREGRESSION_H_

#include "MathMatrix.h"
#include "MathCholesky.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

//use Wald statistics
class LinearRegression{
  public:
    LinearRegression();
    ~LinearRegression();

    bool FitLinearModel(Matrix & X, Vector & y); // return false if not converging
    /// double GetDeviance(Matrix & X, Vector & y);
	Vector & GetAsyPvalue();
	Vector & GetCovEst() {return this->B;} ; // (X'X)^{-1} X'Y
	Matrix & GetCovB() {return this->covB;};
    Vector & GetPredicted() { return this->predict;};
    Vector & GetResiduals() { return this->residuals;};
    double GetSigma2() {return this->sigma2;};
    Vector B;       // coefficient vector
    Matrix covB;    // coefficient covariance matrix
  private:
	Vector pValue;   // 
    Vector residuals;  // Y - X' \hat(beta)
    Vector predict;  // Y - X' \hat(beta)
    Matrix XtXinv;   // (X'X)^ {-1}
    Cholesky chol;
    double sigma2; // \hat{\sigma^2} MLE
};

class LinearRegressionScoreTest{
  public:
    LinearRegressionScoreTest();
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

    double getPvalue() const {return this->pvalue;};

  private:
	void splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol); 
	double pvalue;
    LinearRegression lr;
};


#endif /* _LINEARREGRESSION_H_ */
