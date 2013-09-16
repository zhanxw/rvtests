/////////////////////////////////////////////////////////////////////
// Rewritten by Xiaowei Zhan
//
// Original code is from:
// mach2dat/LogisticRegression.h
// (c) 2008 Yun Li
//
// March 15, 2008
//

#include "LogisticRegression.h"
#include "Eigen/Core"
#include <Eigen/Cholesky>
#include "EigenMatrixInterface.h"
// #include "MathSVD.h"
// #include "MathCholesky.h"
// #include "StringHash.h"
// #include "MathStats.h"

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

class WorkingData {
 public:
  Eigen::MatrixXf X;
  Eigen::VectorXf r;    // residual
  Eigen::VectorXf eta;  // X * beta
  Eigen::VectorXf p;
  Eigen::VectorXf V;    // p * (1-p)
  Eigen::MatrixXf D;  //  X' V X
  Eigen::MatrixXf covB; // (X' V X)^(-1)
  Eigen::VectorXf beta;
  Eigen::VectorXf delta_beta;

  Eigen::VectorXf y;
  Eigen::VectorXf succ;
  Eigen::VectorXf total;
};

LogisticRegression::LogisticRegression()
{
  this->w = new WorkingData;
}

LogisticRegression::~LogisticRegression()
{
  if (this->w) {
    delete this->w;
    this->w = NULL;
  }
}

double LogisticRegression::GetDeviance() {
  double ll = 0.0;
  if (this->w->y.size()) {
    ll = (
        (this->w->y.array() * this->w->p.array().log()) +
        ((1. - this->w->y.array()) * (1.0 - this->w->p.array()).log())
          ).sum();
  } else {
    ll = (
        (this->w->succ.array() * this->w->p.array().log()) +
        ((this->w->total- this->w->succ).array() * (1.0 - this->w->p.array()).log())
          ).sum();
  }

  double deviance = -2.0 * ll;
  return deviance;
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
    double Zstat = B[i] * B[i]/(covB[i][i]);
    // pValue[i] = ndist(Zstat);
    // if (pValue[i] >= 0.5){
    // 	pValue[i] = 2*(1-pValue[i]);
    // } else pValue[i] = 2*pValue[i];
    pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
  }
  return(pValue);
}


void LogisticRegression::Reset(Matrix& X){
  int nr = X.rows;
  int nc = X.cols;

  B.Dimension(nc);
  B.Zero();

  covB.Dimension(nc, nc);
  covB.Zero();

  pValue.Dimension(nc);
  pValue.Zero();

  p.Dimension(nr);
  V.Dimension(nr);
  p.Zero();
  V.Zero();

  this->w->X.setZero(nr, nr);
  this->w->beta.setZero(nc);
  this->w->eta.setZero(nr);
  this->w->p.setZero(nr);
  this->w->V.setZero(nr);
  this->w->D.setZero(nc, nc);
  this->w->covB.setZero(nc, nc);
  this->w->beta.setZero(nc);
  this->w->delta_beta.setZero(nc);

  this->w->y.setZero(0);
  this->w->succ.setZero(0);
  this->w->total.setZero(0);

  // W.Dimension(nr);
  // W.Zero();

  // residuals.Dimension(nr);

  // deltaB.Dimension(nc);

  // D.Dimension(nc, nc);

  // Dinv.Dimension(nc, nc);
  // Dtwo.Dimension(nc, nr);
  // XtV.Dimension(nc, nr);

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

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & succ, Vector& total, int nrrounds) {
  this-> Reset(X);

  G_to_Eigen(X, &this->w->X);
  G_to_Eigen(succ, &this->w->succ);
  G_to_Eigen(total, &this->w->total);

  int rounds = 0;
  double lastDeviance, currentDeviance;
  // Newton-Raphson
  while (rounds < nrrounds) {
    // beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
    this->w->eta = this->w->X * this->w->beta;
    this->w->p = (-this->w->eta.array().exp() + 1.0).inverse();
    this->w->V = this->w->p.array() * (1.0 - this->w->p.array()) * this->w->total.array();

    this->w->D = this->w->X.transpose() * this->w->V.asDiagonal() * this->w->X; // X' V X
    this->w->r = this->w->X.transpose() * (this->w->succ.array() - this->w->total.array() * this->w->p.array()).matrix(); // X' (y-mu)

    this->w->delta_beta = this->w->D.eval().llt().solve(this->w->r);
    if (!(this->w->D * this->w->delta_beta).isApprox(this->w->r, 1e-6)) {
      // cannot inverse
      return false;
    }
    this->w->beta += this->w->delta_beta;
    currentDeviance = this->GetDeviance();
    if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-6)
    { // converged!
      rounds = 0;
      break;
    }
    lastDeviance = currentDeviance;
    rounds ++;
  }
  if (rounds == nrrounds)
  {
    printf("Not enough iterations!");
    return false;
  }
  this->w->covB = this->w->D.eval().llt().solve(Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows()));

  Eigen_to_G(this->w->beta, &B);
  Eigen_to_G(this->w->covB, &covB);
  Eigen_to_G(this->w->p, &p);
  Eigen_to_G(this->w->V, &V);

  return true;

  //   for (int i = 0; i < X.rows; i++)
  //   {
  //     double d = 0;
  //     for (int k = 0; k < X.cols; k ++)
  //     {
  //       d += B[k] * X[i][k]; // \eta,
  //     } // for parameter beta-k

  //     p[i] = 1.0 / (1.0 + exp(-d)); // \mu = E(prob)
  //     V[i] = p[i] * (1 - p[i]);
  //     W[i] = total[i] * V[i]; // weight
  //   } // for observation i

  //   // The first part: solve / inverse
  //   D.Zero();

  //   for (int k = 0; k < X.cols; k++)
  //   {
  //     for (int l = k; l < X.cols; l++)
  //     {
  //       double Dentry = 0.0;
  //       for (int i = 0; i < X.rows; i++)
  //         Dentry += X[i][k] * W[i] * X[i][l];
  //       D[k][l] = D[l][k] = Dentry;
  //     }
  //   }

  //   Dinv.Zero();
  //   Dinv = D;
  //   // SVD svd;
  //   // svd.SVDInverse(Dinv);
  //   if (!chol.TryDecompose(D))
  //     return false;
  //   chol.Decompose(D);
  //   chol.Invert();
  //   Dinv = chol.inv; // (X' W X)^{-1}

  //   // Accumulate up to the second part: multiply by t(X)
  //   Dtwo.Zero();

  //   for (int k = 0; k < X.cols; k++)
  //     for (int i = 0; i < X.rows; i++)
  //       for (int l = 0; l < X.cols; l++)
  //         Dtwo[k][i] += Dinv[k][l] * X[i][l];

  //   // The last part, residuals
  //   residuals.Zero();
  //   for (int i = 0; i < X.rows; i++)
  //     residuals[i] = succ[i] - total[i] * p[i];

  //   deltaB.Zero();
  //   for (int k = 0; k < X.cols; k++)
  //     for (int i = 0; i < X.rows; i++)
  //       deltaB[k] += Dtwo[k][i] * residuals[i];

  //   // update beta's and check for convergence
  //   for (int k = 0; k < X.cols; k++)
  //   {
  //     B[k] += deltaB[k];
  //   }

  //   //printf ("%d %-11.4g\n", rounds, delta);
  //   currentDeviance = this->GetDeviance(X,succ, total);
  //   if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-6)
  //   {
  //     rounds = 0;
  //     break;
  //   }
  //   lastDeviance = currentDeviance;
  //   rounds ++;
  // } // Newton-Raphson iterations

  // if (rounds == nrrounds)
  //   return false;


  // // obtain covariance matrix to perform Wald test
  // // covB = solve(t(X)%*%V%*%X)

  // // Transpose X and multiply by diagonal V
  // //XtV.Zero();
  // for (int i = 0; i < X.rows; i++)
  //   for (int k = 0; k < X.cols; k++)
  //     XtV[k][i] = X[i][k] * W[i];

  // //  use this is simpler one line to replace following:
  // // 	Matrix XtVX;
  // // 	XtVX.Dimension(X.cols, X.cols);
  // // 	XtVX.Zero();
  // // 	XtVX.Product(XtV,X);

  // // 	covB.Zero();
  // // 	covB = XtVX;
  // covB.Product(XtV, X);

  // //SVD svd;
  // //svd.SVDInverse(covB);
  // //Cholesky chol;
  // if (!chol.TryDecompose(covB))
  //   return false;
  // chol.Decompose(covB);
  // chol.Invert();
  // covB = chol.inv;

  // return true;
}

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & y, int nrrounds)
{
  this-> Reset(X);

  G_to_Eigen(X, &this->w->X);
  G_to_Eigen(y, &this->w->y);

  int rounds = 0;
  double lastDeviance, currentDeviance;
  while (rounds < nrrounds) {
    this->w->eta = this->w->X * this->w->beta;
    this->w->p = (1.0 + (-this->w->eta.array()).exp()).inverse();
    this->w->V = this->w->p.array() * (1.0 - this->w->p.array());
    this->w->D = this->w->X.transpose() * this->w->V.asDiagonal() * this->w->X; // X' V X
    this->w->r = this->w->X.transpose() * (this->w->y - this->w->p); // X' (y-mu)

    this->w->delta_beta = this->w->D.eval().llt().solve(this->w->r);
    if (!(this->w->D * this->w->delta_beta).isApprox(this->w->r, 1e-6)) {
      // cannot inverse
      return false;
    }
    this->w->beta += this->w->delta_beta;
    currentDeviance = this->GetDeviance();
    if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-6)
    { // converged!
      rounds = 0;
      break;
    }
    lastDeviance = currentDeviance;
    rounds ++;
  }
  if (rounds == nrrounds)
  {
    printf("Not enough iterations!");
    return false;
  }
  this->w->covB = this->w->D.eval().llt().solve(Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows()));

  Eigen_to_G(this->w->beta, &B);
  Eigen_to_G(this->w->covB, &covB);
  Eigen_to_G(this->w->p, &p);
  Eigen_to_G(this->w->V, &V);

  return true;


#if 0
  {
    this-> Reset(X);
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

  }
  return true;
#endif
}
