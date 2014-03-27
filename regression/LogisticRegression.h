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
/* #if 1 */
/* #include "MathCholesky.h" */
/* #endif */
/* #include "StringHash.h" */
/* #include "StringArray.h" */
/* #include <cmath> */
/* #include "MathStats.h" */
/* #include "MathSVD.h" */


class WorkingData; // store temporary data structure

//use Wald statistics
class LogisticRegression
{
 public:
  LogisticRegression();
  ~LogisticRegression();

  // main function
  // @return true if success; false if not converging
  bool FitLogisticModel(Matrix & X, Matrix & y, int rnrounds); 
  bool FitLogisticModel(Matrix & X, Vector & y, int rnrounds); 
  bool FitLogisticModel(Matrix & X, Vector & succ, Vector& total, int nrrounds);

  // alias simplified functions
  bool Fit(Matrix & X, Matrix & y) {
    return this->FitLogisticModel(X, y, 100);
  }
  bool Fit(Matrix & X, Vector & y) {
    return this->FitLogisticModel(X, y, 100);
  }
  bool Fit(Matrix & X, Vector & succ, Vector& total) {
    return this->FitLogisticModel(X, succ, total, 100);
  }

  // obtain results
  double GetDeviance(Matrix & X, Vector & y);
  double GetDeviance(Matrix & X, Vector & succ, Vector& total);
  Vector& GetAsyPvalue();
  Vector& GetCovEst()    {return this->B;}; // coef estimation of the model
  Vector& GetPredicted() {return this->p;}; // predicted probability \hat{p}
  Vector& GetVariance()  {return this->V;}; // predicted variance ( \hat{p} * (1- \hat{p})
  Matrix& GetCovB()      {return this->covB;} ;

  //todo: Fitxxx() will call Reset(), and thus make this not useful...
  //      may change Reset() function
  void SetInitialCovEst(Vector& initB) { this->B = initB;} ; // set initial value of B, that may speed estimation up if this initial value is close to estimated results.
  void Reset(Matrix& X); // get everything cleared

 private:
  double GetDeviance();

 private:
  // dont' copy.
  LogisticRegression(const LogisticRegression& l);
  LogisticRegression& operator=(const LogisticRegression& l);

 private:
  Vector B;             // coefficient vector
  Matrix covB;          // coefficient covariance matrix
  Vector pValue;        // pvalues
  Vector p;             // p: estimted prob;
  Vector V;             // V: p(1-p) ;
#if 0
  Vector W;             // W n*p*(1-p)
  Vector residuals;
  Vector deltaB;
  Matrix D;
  // Matrix testSummary;
  Matrix Dinv;
  Cholesky chol;
  Matrix Dtwo;
  Matrix XtV;
#endif
  WorkingData* w;       // holding temporary caluclation results
};
#endif
