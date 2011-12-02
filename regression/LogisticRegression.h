//////////////////////////////////////////////////////////////////////
// Adopted and improved by Xiaowei Zhan, Youna Hu
// Score Test added
//
// mach2dat/LogisticRegression.h
// (c) 2008 Yun Li
//
// March 15, 2008
//

#ifndef __LOGISTIC_REGRESSION_H__
#define __LOGISTIC_REGRESSION_H__

#include "MathMatrix.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "StringArray.h"
#include <cmath>
#include "MathStats.h"
#include "MathSVD.h"

//use Wald statistics
class LogisticRegression
{
  public:
    LogisticRegression();
    ~LogisticRegression();

    bool FitLogisticModel(Matrix & X, Vector & y, int rnrounds); // return false if not converging
    bool FitLogisticModel(Matrix & X, Vector & succ, Vector& total, int nrrounds);
    double GetDeviance(Matrix & X, Vector & y);
    double GetDeviance(Matrix & X, Vector & succ, Vector& total);
	Vector & GetAsyPvalue();
	Vector & GetCovEst();
	Matrix & GetCovB();
    void reset(Matrix& X); // get everything cleared

    Vector B;       // coefficient vector
    Matrix covB;    // coefficient covariance matrix

  private:
	Vector pValue;
	Vector p, V, W;
    Vector residuals;
    Vector deltaB;
    Matrix D;
	Matrix testSummary;
    Matrix Dinv;
    Cholesky chol;
    Matrix Dtwo;
	Matrix XtV;

};

class LogisticRegressionScoreTest{
  public:
    LogisticRegressionScoreTest();
	bool FitLogisticModel(Matrix &X, Vector &y, int colToTest, int nRound);
    bool FitNullModel(Matrix& Xnull, Vector& y, int nRound);
    bool TestCovariate(Matrix& X, int colToTest);
    double getPvalue() const {return this->pvalue;};
  private:
	void splitMatrix(Matrix& x, int col, Matrix& xnull, Vector& xcol); 
	double pvalue;
};
#endif



