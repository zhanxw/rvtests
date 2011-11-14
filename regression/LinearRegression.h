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

    Vector B;       // coefficient vector
    Matrix covB;    // coefficient covariance matrix
  private:
	Vector pValue;   // 
    Vector residuals;  // Y - X' \hat(beta)
    Vector predict;  // Y - X' \hat(beta)
    Matrix XtXinv;   // (X'X)^ {-1}
    Cholesky chol;
};


#endif /* _LINEARREGRESSION_H_ */
