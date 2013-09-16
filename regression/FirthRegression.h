#ifndef __FIRTH_REGRESSION_H__
#define __FIRTH_REGRESSION_H__

#include "MathMatrix.h"
class WorkingData; // store temporary data structure

//use Wald statistics
class FirthRegression
{
 public:
  FirthRegression();
  ~FirthRegression();

  // main function
  bool FitFirthModel(Matrix & X, Matrix & y, int rnrounds); // return false if not converging
  bool FitFirthModel(Matrix & X, Vector & y, int rnrounds); // return false if not converging
  bool FitFirthModel(Matrix & X, Vector & succ, Vector& total, int nrrounds);

  // alias simplified functions
  bool Fit(Matrix & X, Matrix & y) {
    return this->FitFirthModel(X, y, 100);
  }
  bool Fit(Matrix & X, Vector & y) {
    return this->FitFirthModel(X, y, 100);
  }
  bool Fit(Matrix & X, Vector & succ, Vector& total) {
    return this->FitFirthModel(X, succ, total, 100);
  }

  // obtain results
  /* double GetDeviance(Matrix & X, Vector & y); */
  /* double GetDeviance(Matrix & X, Vector & succ, Vector& total); */
  Vector & GetAsyPvalue();
  Vector & GetCovEst()    {return this->B;}; // coef estimation of the model
  Vector & GetPredicted() {return this->p;}; // predicted probability \hat{p}
  Vector & GetVariance()  {return this->V;}; // predicted variance ( \hat{p} * (1- \hat{p})
  Matrix & GetCovB()      {return this->covB;} ;

  void SetInitialCovEst(Vector& initB) { this->B = initB;} ; // set initial value of B, that may speed estimation up if this initial value is close to estimated results.
  void Reset(Matrix& X); // get everything cleared

  /* private: */
  /*   double GetDeviance(); */
 private:
  // dont' copy.
  FirthRegression(const FirthRegression& l);
  FirthRegression& operator=(const FirthRegression& l);

 private:
  Vector B;             // coefficient vector
  Matrix covB;          // coefficient covariance matrix
  Vector pValue;        // pvalues
  Vector p;             // p: estimted prob;
  Vector V;             // V: p(1-p) ;
  WorkingData* w;       // holding temporary caluclation results
};
#endif
