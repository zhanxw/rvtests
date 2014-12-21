#include "LinearRegression.h"
#include "gsl/gsl_cdf.h"

bool LinearRegression::FitLinearModel(Matrix & X, Matrix & y) {
  if (y.cols != 1) {
    fprintf(stderr, "%s:%d Use first column of y\n", __FILE__, __LINE__);
  }
  Vector v(y.rows);
  for (int i = 0; i < y.rows; ++i){
    v[i] = y[i][0];
  }
  return this->FitLinearModel(X, v);
}

bool LinearRegression::FitLinearModel(Matrix & X, Vector & y){
  Matrix Xt;
  Xt.Transpose(X);

  Matrix XtX;
  XtX.Product(Xt, X);
  if (!this->chol.TryDecompose(XtX))
    return false;
  chol.Decompose(XtX);
  chol.Invert();
  this->XtXinv = chol.inv;

  Vector tmp = y;
  tmp.Product(Xt, y);
  this->B.Product(this->XtXinv, tmp); // beta = (XtX)^{-1} Xt Y

  this->predict.Product(X, this->B);
  this->residuals = y;
  this->residuals.Subtract(this->predict);

  this->sigma2 = 0.0;
  for (int i = 0; i < this->residuals.Length(); i++){
    sigma2 += (this->residuals[i]) * (this->residuals[i]);
  }
  sigma2 /= y.Length(); // MLE estimates of sigma2

  this->covB = this->XtXinv;
  this->covB.Multiply(sigma2);
  return true;
};

Vector& LinearRegression::GetAsyPvalue(){
  int numCov = B.Length();
  pValue.Dimension(B.Length());
  for (int i = 0; i < numCov; i ++){
    double Zstat = B[i]/sqrt(covB[i][i]);
    // pValue[i] = ndist(Zstat);
    // if (pValue[i] >= 0.5){
    //      pValue[i] = 2*(1-pValue[i]);
    // } else pValue[i] = 2*pValue[i];
    Zstat *= Zstat;
    pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
  }
  return(pValue);
}

bool LinearRegression::calculateResidualMatrix(Matrix& X, Matrix* out) {
  if (!calculateHatMatrix(X, out)) return false;
  Matrix& m = *out;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      if (i == j) {
        m[i][j] = 1.0 - m[i][j];
      } else {
        m[i][j] = - m[i][j];
      }
    }
  }
  return true;
}
bool LinearRegression::calculateHatMatrix(Matrix& X, Matrix* out) {
  Matrix Xt;
  Xt.Transpose(X);

  Matrix XtX;
  XtX.Product(Xt, X);
  if (!this->chol.TryDecompose(XtX))
    return false;
  chol.Decompose(XtX);
  chol.Invert();
  this->XtXinv = chol.inv;

  Matrix tmp;
  tmp.Product(XtXinv, Xt);
  (*out).Product(X, tmp);
  return true;
}
