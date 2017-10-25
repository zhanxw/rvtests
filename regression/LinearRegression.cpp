#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include "LinearRegression.h"

#include <stdio.h>  // printf
#include "third/eigen/Eigen/Cholesky"
#include "third/gsl/include/gsl/gsl_cdf.h"

bool LinearRegression::FitLinearModel(Matrix& X, Matrix& y) {
  if (y.cols != 1) {
    fprintf(stderr, "%s:%d Use first column of y\n", __FILE__, __LINE__);
  }
  Vector v(y.rows);
  for (int i = 0; i < y.rows; ++i) {
    v[i] = y(i, 0);
  }
  return this->FitLinearModel(X, v);
}

bool LinearRegression::FitLinearModel(Matrix& X, Vector& y) {
  // Matrix Xt;
  // Xt.Transpose(X);

  // Matrix XtX;
  // XtX.Product(Xt, X);
  // if (!this->chol.TryDecompose(XtX)) return false;
  // chol.Decompose(XtX);
  // chol.Invert();
  // this->XtXinv = chol.inv;
  XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  DECLARE_EIGEN_MATRIX(XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e)
                 .llt()
                 .solve(Eigen::MatrixXd::Identity(X_e.cols(), X_e.cols()));

  // Vector tmp = y;
  // tmp.Product(Xt, y);
  // this->B.Product(this->XtXinv, tmp);  // beta = (XtX)^{-1} Xt Y
  B.Dimension(X.cols, 1);
  DECLARE_EIGEN_VECTOR(B, B_e);
  DECLARE_EIGEN_CONST_VECTOR(y, y_e);
  B_e = XtXinv_e * X_e.transpose() * y_e;

  // this->predict.Product(X, this->B);
  // this->residuals = y;
  // this->residuals.Subtract(this->predict);
  this->predict.Dimension(X.rows, 1);
  this->residuals.Dimension(X.rows, 1);
  DECLARE_EIGEN_VECTOR(this->predict, predict_e);
  DECLARE_EIGEN_VECTOR(this->residuals, resid_e);
  predict_e = X_e * B_e;
  resid_e = y_e - predict_e;

  // this->sigma2 = 0.0;
  // for (int i = 0; i < this->residuals.Length(); i++) {
  //   sigma2 += (this->residuals[i]) * (this->residuals[i]);
  // }
  // sigma2 /= y.Length();  // MLE estimates of sigma2
  this->sigma2 = resid_e.squaredNorm() / y_e.size();

  // this->covB = this->XtXinv;
  // this->covB.Multiply(sigma2);
  this->covB.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(this->covB, covB_e);
  covB_e = XtXinv_e * sigma2;

  return true;
};

Vector& LinearRegression::GetAsyPvalue() {
  int numCov = B.Length();
  pValue.Dimension(B.Length());
  for (int i = 0; i < numCov; i++) {
    double Zstat = B[i] / sqrt(covB(i, i));
    // pValue[i] = ndist(Zstat);
    // if (pValue[i] >= 0.5){
    //      pValue[i] = 2*(1-pValue[i]);
    // } else pValue[i] = 2*pValue[i];
    Zstat *= Zstat;
    pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
  }
  return (pValue);
}

bool LinearRegression::calculateResidualMatrix(Matrix& X, Matrix* out) {
  if (!calculateHatMatrix(X, out)) return false;
  Matrix& m = *out;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      if (i == j) {
        m(i, j) = 1.0 - m(i, j);
      } else {
        m(i, j) = -m(i, j);
      }
    }
  }
  return true;
}
bool LinearRegression::calculateHatMatrix(Matrix& X, Matrix* out) {
  // Matrix Xt;
  // Xt.Transpose(X);

  // Matrix XtX;
  // XtX.Product(Xt, X);
  // if (!this->chol.TryDecompose(XtX)) return false;
  // chol.Decompose(XtX);
  // chol.Invert();
  // this->XtXinv = chol.inv;
  DECLARE_EIGEN_CONST_MATRIX(X, X_e);
  this->XtXinv.Dimension(X.cols, X.cols);
  DECLARE_EIGEN_MATRIX(this->XtXinv, XtXinv_e);
  XtXinv_e = (X_e.transpose() * X_e)
                 .llt()
                 .solve(Eigen::MatrixXd::Identity(X_e.cols(), X_e.cols()));

  // Matrix tmp;
  // tmp.Product(XtXinv, Xt);
  // (*out).Product(X, tmp);

  DECLARE_EIGEN_MATRIX((*out), out_e);
  out_e = X_e * XtXinv_e * X_e.transpose();
  return true;
}
