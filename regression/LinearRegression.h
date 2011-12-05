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
    Vector B;       // coefficient vector
    Matrix covB;    // coefficient covariance matrix
  private:
	Vector pValue;   // 
    Vector residuals;  // Y - X' \hat(beta)
    Vector predict;  // Y - X' \hat(beta)
    Matrix XtXinv;   // (X'X)^ {-1}
    Cholesky chol;
};

class LinearRegressionScoreTest{
  public:
    LinearRegressionScoreTest();
	bool FitLinearModel(Matrix &X, Vector &y, int colToTest);

    bool FitNullModel(Matrix& Xnull, Vector& y);
    bool TestCovariate(Matrix& Xnull, Vector& y, Vector& Xcol);
    
    // fit y~1+ beta*x  (no covariate)
    bool TestCovariate(Vector& x, Vector& y);

    double getPvalue() const {return this->pvalue;};

  private:
	void splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol); 
	double pvalue;
    LinearRegression lr;
};


#endif /* _LINEARREGRESSION_H_ */
