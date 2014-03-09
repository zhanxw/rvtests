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
#include <math.h>
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
  Eigen::MatrixXd X;
  Eigen::VectorXd r;    // residual
  Eigen::VectorXd eta;  // X * beta
  Eigen::VectorXd p;
  Eigen::VectorXd V;    // p * (1-p)
  Eigen::MatrixXd D;  //  X' V X
  Eigen::MatrixXd covB; // (X' V X)^(-1)
  Eigen::VectorXd beta;
  Eigen::VectorXd delta_beta;

  Eigen::VectorXd y;
  Eigen::VectorXd succ;
  Eigen::VectorXd total;
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
  // int i1x, i1y, i2x, i2y;
  // fprintf(stderr, "min[%d][%d] = %g, max[%d][%d] = %g\n",
  //         i1x, i1y, w->p.minCoeff(&i1x, &i1y),
  //         i2x, i2y, w->p.maxCoeff(&i2x, &i2y));
  
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
  const int nSample = this->w->X.rows();
  // Newton-Raphson
  while (rounds < nrrounds) {
    // beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
    this->w->eta = this->w->X * this->w->beta;
    this->w->p = ((-this->w->eta.array()).exp() + 1.0).inverse();
    this->w->V = this->w->p.array() * (1.0 - this->w->p.array()) * this->w->total.array();

    this->w->D = this->w->X.transpose() * this->w->V.asDiagonal() * this->w->X; // X' V X
    this->w->r = this->w->X.transpose() * (this->w->succ.array() - this->w->total.array() * this->w->p.array()).matrix(); // X' (y-mu)

    this->w->delta_beta = this->w->D.eval().llt().solve(this->w->r);
    // const double rel = (this->w->D * this->w->delta_beta - this->w->r).norm() / this->w->r.norm();
    // if ( this->w->r.norm() >0 && rel > 1e-6) {
    double norm = (this->w->D * this->w->delta_beta - this->w->r).norm();
    if (norm  > 1e-3 * nSample) {
      // cannot inverse
      return false;
    }
    this->w->beta += this->w->delta_beta;
    currentDeviance = this->GetDeviance();
    if (rounds >1 && fabs(currentDeviance - lastDeviance) < 1e-3)
    { // converged!
      rounds = 0;
      break;
    }
    if (std::fpclassify(currentDeviance) != FP_NORMAL) {
      // probably separation happens
      return false;
    }
    lastDeviance = currentDeviance;
    rounds ++;
  }
  if (rounds == nrrounds)
  {
    printf("[ %s:%d ] Not enough iterations!", __FILE__, __LINE__);
    return false;
  }
  this->w->covB = this->w->D.eval().llt().solve(Eigen::MatrixXd::Identity(this->w->D.rows(), this->w->D.rows()));

  Eigen_to_G(this->w->beta, &B);
  Eigen_to_G(this->w->covB, &covB);
  Eigen_to_G(this->w->p, &p);
  Eigen_to_G(this->w->V, &V);

  return true;
}

bool LogisticRegression::FitLogisticModel(Matrix & X, Vector & y, int nrrounds)
{
  this-> Reset(X);

  G_to_Eigen(X, &this->w->X);
  G_to_Eigen(y, &this->w->y);

  int rounds = 0;
  double lastDeviance, currentDeviance;
  // const int nSample = this->w->X.rows();
  while (rounds < nrrounds) {
    this->w->eta = this->w->X * this->w->beta;
    this->w->p = (1.0 + (-this->w->eta.array()).exp()).inverse();
    this->w->V = this->w->p.array() * (1.0 - this->w->p.array());
    this->w->D = this->w->X.transpose() * this->w->V.asDiagonal() * this->w->X; // X' V X
    this->w->r = this->w->X.transpose() * (this->w->y - this->w->p); // X' (y-mu)

    this->w->delta_beta = this->w->D.eval().llt().solve(this->w->r);

    // // output different norms, see which works
    // fprintf(stderr, "norm[1] = %g\n", (this->w->D * this->w->delta_beta - this->w->r).norm());
    // fprintf(stderr, "norm[2] = %s\n", (this->w->D * this->w->delta_beta).isApprox(this->w->r, 1e-3) ? "true": "false");
    // fprintf(stderr, "norm[2] = %s\n", (this->w->D * this->w->delta_beta).isApprox(this->w->r, 1e-6) ? "true": "false");
    
    // const double rel = (this->w->D * this->w->delta_beta - this->w->r).norm() / this->w->r.norm();
    // fprintf(stderr, "rel = %g, norm1 = %g, norm2 = %g\n", rel, (this->w->D * this->w->delta_beta - this->w->r).norm(), this->w->r.norm());

    // if ( this->w->r.norm() >0 && rel > 1e-6) {
    //   // cannot inverse
    //   // return false;
    // }
#if 0
    // D = X' V X (positive definite), so cholesky decomposition should always work
    double norm = (this->w->D * this->w->delta_beta - this->w->r).norm();
    if (norm > 1e-3 * nSample) { // very practical choice...
      // cannot inverse
      return false;
    }
#endif 
    this->w->beta += this->w->delta_beta;
    currentDeviance = this->GetDeviance();

    // printf("%lf, %lf, %lf, %s\n", currentDeviance, lastDeviance, std::fabs(currentDeviance - lastDeviance),
    //        (rounds >1 && std::fabs(currentDeviance - lastDeviance) < 1e-3)?"true":"false");
    
    if (rounds >1 && std::fabs(currentDeviance - lastDeviance) < 1e-3)
    { // converged!
      rounds = 0;
      break;
    }
    if (std::fpclassify(currentDeviance) != FP_NORMAL) {
      // probably separation happens
      return false;
    }
    lastDeviance = currentDeviance;
    rounds ++;
  }
  if (rounds == nrrounds)
  {
    printf("[ %s:%d ] Not enough iterations!\n", __FILE__, __LINE__);
    printf("%lf, %lf, %lf, %s\n", currentDeviance, lastDeviance, std::fabs(currentDeviance - lastDeviance),
           (rounds >1 && std::fabs(currentDeviance - lastDeviance) < 1e-3)?"true":"false");
    return false;
  }
  this->w->covB = this->w->D.eval().llt().solve(Eigen::MatrixXd::Identity(this->w->D.rows(), this->w->D.rows()));

  Eigen_to_G(this->w->beta, &B);
  Eigen_to_G(this->w->covB, &covB);
  Eigen_to_G(this->w->p, &p);
  Eigen_to_G(this->w->V, &V);

  return true;
}
