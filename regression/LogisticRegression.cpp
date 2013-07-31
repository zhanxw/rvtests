/////////////////////////////////////////////////////////////////////
// Adopted by Xiaowei Zhan
//
// mach2dat/LogisticRegression.h
// (c) 2008 Yun Li
//
// March 15, 2008
//

#include "LogisticRegression.h"
#include "MathSVD.h"
#include "MathCholesky.h"
#include "StringHash.h"
#include "MathStats.h"

#include "gsl/gsl_cdf.h"

#ifdef DEBUG
// for debug usage
void printToFile(Vector& v, String fn, int index) {
  String n;
  n.printf("%s%d",fn.c_str(), index);
  FILE* fp = fopen(n.c_str(), "wt");
  for (int i = 0; i< v.Length(); i++)
    fprintf(fp, "%lf\n", v[i]);
  fclose(fp);
}
// for debug usage
void printToFile(Matrix& m, String fn, int index) {
  String n;
  n.printf("%s%d",fn.c_str(), index);
  FILE* fp = fopen(n.c_str(), "wt");
  for (int i = 0; i<m.rows; i++) {
    for (int j = 0; j< m.cols; j++) {
      fprintf(fp, "%lf\t", m[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
#endif

LogisticRegression::LogisticRegression()
{
}

LogisticRegression::~LogisticRegression()
{
}

double LogisticRegression::GetDeviance(Matrix & X, Vector & y)
{
  double ll = 0.0;
  for (int i = 0; i < X.rows; i++)
  {
    double t = 0.0;
    for (int j = 0; j < X.cols; j++)
      t += B[j] * X[i][j];
    double yhat = 1 / (1 + exp(-t));

    ll += y[i] == 1 ? log (yhat) : log(1 - yhat);
  }

  double deviance = -2.0 * ll;
  return deviance;
}

double LogisticRegression::GetDeviance(Matrix & X, Vector & succ, Vector& total)
{
  double ll = 0.0;
  for (int i = 0; i < X.rows; i++)
  {
    double t = 0.0;
    for (int j = 0; j < X.cols; j++)
      t += B[j] * X[i][j];
    double yhat = 1 / (1 + exp(-t));

    ll += succ[i] * log(yhat) + (total[i] - succ[i]) * log(1-yhat);
  }

  double deviance = -2.0 * ll;
  return deviance;
}

Vector & LogisticRegression::GetAsyPvalue(){
  int numCov = B.Length();
  pValue.Dimension(B.Length());
  for (int i = 0; i < numCov; i ++){
    double Zstat = B[i]/sqrt(covB[i][i]);
    // pValue[i] = ndist(Zstat);
    // if (pValue[i] >= 0.5){
    // 	pValue[i] = 2*(1-pValue[i]);
    // } else pValue[i] = 2*pValue[i];
    Zstat *= Zstat;
    pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
  }
  return(pValue);
}


void LogisticRegression::reset(Matrix& X){
  B.Dimension(X.cols);
  B.Zero();

  covB.Dimension(X.cols, X.cols);
  covB.Zero();

  p.Dimension(X.rows);
  V.Dimension(X.rows);
  W.Dimension(X.rows);
  p.Zero();
  V.Zero();
  W.Zero();

  residuals.Dimension(X.rows);

  deltaB.Dimension(X.cols);

  D.Dimension(X.cols, X.cols);

  Dinv.Dimension(X.cols, X.cols);
  Dtwo.Dimension(X.cols, X.rows);
  XtV.Dimension(X.cols, X.rows);

}

bool LogisticRegression::FitLogisticModel(Matrix & X, Matrix & y, int rnrounds) {
  if (y.cols != 1) {
    fprintf(stderr, "%s:%d Use first column of y\n", __FILE__, __LINE__);
  }
  Vector v(X.rows);
  for (int i = 0; i < X.rows; ++i){
    v[i] = y[i][0];
  }
  return this->FitLogisticModel(X, v, rnrounds);
};

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & succ, Vector& total, int nrrounds){
  this-> reset(X);
  int rounds = 0;
  double lastDeviance, currentDeviance;
  // Newton-Raphson
  while (rounds < nrrounds)
    //while (true)
  {
    // beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
    for (int i = 0; i < X.rows; i++)
    {
      double d = 0;
      for (int k = 0; k < X.cols; k ++)
      {
        d += B[k] * X[i][k]; // \eta,
      } // for parameter beta-k

      p[i] = 1.0 / (1.0 + exp(-d)); // \mu = E(prob)
      V[i] = p[i] * (1 - p[i]);
      W[i] = total[i] * V[i]; // weight
    } // for observation i

    // The first part: solve / inverse
    D.Zero();

    for (int k = 0; k < X.cols; k++)
    {
      for (int l = k; l < X.cols; l++)
      {
        double Dentry = 0.0;
        for (int i = 0; i < X.rows; i++)
          Dentry += X[i][k] * W[i] * X[i][l];
        D[k][l] = D[l][k] = Dentry;
      }
    }

    Dinv.Zero();
    Dinv = D;
    // SVD svd;
    // svd.SVDInverse(Dinv);
    if (!chol.TryDecompose(D))
      return false;
    chol.Decompose(D);
    chol.Invert();
    Dinv = chol.inv; // (X' W X)^{-1}

    // Accumulate up to the second part: multiply by t(X)
    Dtwo.Zero();

    for (int k = 0; k < X.cols; k++)
      for (int i = 0; i < X.rows; i++)
        for (int l = 0; l < X.cols; l++)
          Dtwo[k][i] += Dinv[k][l] * X[i][l];

    // The last part, residuals
    residuals.Zero();
    for (int i = 0; i < X.rows; i++)
      residuals[i] = succ[i] - total[i] * p[i];

    deltaB.Zero();
    for (int k = 0; k < X.cols; k++)
      for (int i = 0; i < X.rows; i++)
        deltaB[k] += Dtwo[k][i] * residuals[i];

    // update beta's and check for convergence
    for (int k = 0; k < X.cols; k++)
    {
      B[k] += deltaB[k];
    }

    //printf ("%d %-11.4g\n", rounds, delta);
    currentDeviance = this->GetDeviance(X,succ, total);
    if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-6)
    {
      rounds = 0;
      break;
    }
    lastDeviance = currentDeviance;
    rounds ++;
  } // Newton-Raphson iterations

  if (rounds == nrrounds)
    return false;


  // obtain covariance matrix to perform Wald test
  // covB = solve(t(X)%*%V%*%X)

  // Transpose X and multiply by diagonal V
  //XtV.Zero();
  for (int i = 0; i < X.rows; i++)
    for (int k = 0; k < X.cols; k++)
      XtV[k][i] = X[i][k] * W[i];

  //  use this is simpler one line to replace following:
  // 	Matrix XtVX;
  // 	XtVX.Dimension(X.cols, X.cols);
  // 	XtVX.Zero();
  // 	XtVX.Product(XtV,X);

  // 	covB.Zero();
  // 	covB = XtVX;
  covB.Product(XtV, X);

  //SVD svd;
  //svd.SVDInverse(covB);
  //Cholesky chol;
  if (!chol.TryDecompose(covB))
    return false;
  chol.Decompose(covB);
  chol.Invert();
  covB = chol.inv;

  return true;
}

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & y, int nrrounds)
{
  this-> reset(X);
  int rounds = 0;
  double lastDeviance, currentDeviance;
  /*double dChecksum;
    for (int b = 0; b < X.cols; b++) {
    dChecksum = 0.0;
    for (int a = 0; a < X.rows; a++) {
    dChecksum += X[a][b];
    }
    printf("%1.15e\n", dChecksum);
    } */

  // Newton-Raphson
  while (rounds < nrrounds)
    //while (true)
  {
    // beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
    for (int i = 0; i < X.rows; i++)
    {
      double d = 0;
      for (int k = 0; k < X.cols; k ++)
      {
        d += B[k] * X[i][k]; // \eta,
      } // for parameter beta-k

      p[i] = 1.0 / (1.0 + exp(-d)); // \mu = E(prob)
      V[i] = p[i] * (1 - p[i]); // weight
    } // for observation i

    // The first part: solve / inverse
    D.Zero();

    for (int k = 0; k < X.cols; k++)
    {
      for (int l = k; l < X.cols; l++)
      {
        double Dentry = 0.0;
        for (int i = 0; i < X.rows; i++)
          Dentry += X[i][k] * V[i] * X[i][l];
        D[k][l] = D[l][k] = Dentry;
      }
    }
    Dinv.Zero();
    Dinv = D;
    // SVD svd;
    // svd.SVDInverse(Dinv);
    //printf("The size of D = (%d,%d)\n",D.rows,D.cols);

    if (!chol.TryDecompose(D)){
#ifndef NDEBUG
      printf("Cannot decompose D!");
      printf("Matrix D \n");
      for (int i = 0; i < D.rows; i ++) {
        for (int j = 0; j < D.cols; j ++) {
          printf("%g\t", D[i][j]);
        }
        printf("\n");
      }
      // printf("Matrix X \n");
      // for (int i = 0; i < X.rows; i ++) {
      //   for (int j = 0; j < X.cols; j ++) {
      //     printf("%g\t", X[i][j]);
      //   }
      //   printf("\n");
      // }
#endif
      return false;
    }
    chol.Decompose(D);
    chol.Invert();
    Dinv = chol.inv; // (X' W X)^{-1}
    // Accumulate up to the second part: multiply by t(X)
    // now Dtwo = (X' W X)^{-1} X'
    Dtwo.Zero();

    for (int k = 0; k < X.cols; k++)
      for (int i = 0; i < X.rows; i++)
        for (int l = 0; l < X.cols; l++)
          Dtwo[k][i] += Dinv[k][l] * X[i][l];

    // The last part, residuals
    residuals.Zero();
    for (int i = 0; i < X.rows; i++)
      residuals[i] = y[i] - p[i];
    // deltaB = (X' W X)^{-1} X' (y-p)
    deltaB.Zero();
    for (int k = 0; k < X.cols; k++)
      for (int i = 0; i < X.rows; i++)
        deltaB[k] += Dtwo[k][i] * residuals[i];

    // update beta's and check for convergence
    for (int k = 0; k < X.cols; k++)
    {
      B[k] += deltaB[k];
    }
    currentDeviance = this->GetDeviance(X,y);
    if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-6)
    {
      rounds = 0;
      break;
    }
    lastDeviance = currentDeviance;
    rounds ++;
  } // Newton-Raphson iterations

  if (rounds == nrrounds)
  {
    printf("Not enough iterations!");
    return false;
  }



  // obtain covariance matrix to perform Wald test
  // covB = solve(t(X)%*%V%*%X)

  // Transpose X and multiply by diagonal V
  //XtV.Zero();
  for (int i = 0; i < X.rows; i++)
    for (int k = 0; k < X.cols; k++)
      XtV[k][i] = X[i][k] * V[i];

  //  use this is simpler one line to replace following:
  // 	Matrix XtVX;
  // 	XtVX.Dimension(X.cols, X.cols);
  // 	XtVX.Zero();
  // 	XtVX.Product(XtV,X);

  // 	covB.Zero();
  // 	covB = XtVX;
  covB.Product(XtV, X);

  //SVD svd;
  //svd.SVDInverse(covB);
  //Cholesky chol;
  if (!chol.TryDecompose(covB)){
    printf("cannot decompose covB");
    return false;
  }

  chol.Decompose(covB);
  chol.Invert();
  covB = chol.inv;

  return true;
}
