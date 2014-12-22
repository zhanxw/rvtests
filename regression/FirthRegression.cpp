#include "FirthRegression.h"

#include "Eigen/Core"
#include <Eigen/Cholesky>
#include "EigenMatrixInterface.h"

#include "gsl/gsl_cdf.h"

#ifndef NDEBUG
#include <iostream>
#include <fstream>
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
void printToFile(Eigen::MatrixXf& m, const char* fn, const char* label) {
  std::ofstream ofs (fn, std::ofstream::out);
  ofs << "#[ " << label << " ]\n";
  ofs << m << "\n";
  ofs.close();
}
void printToFile(Eigen::VectorXf& v, const char* fn, const char* label) {
  std::ofstream ofs (fn, std::ofstream::out);
  ofs << "#[ " << label << " ]\n";
  ofs << v << "\n";
  ofs.close();
}
#endif

class WorkingData {
 public:
  Eigen::MatrixXf X;
  Eigen::VectorXf r;    // residual
  Eigen::VectorXf eta;  // X * beta
  Eigen::VectorXf p;    // p = 1/(1+exp(-eta))
  Eigen::VectorXf V;    // p * (1-p)
  Eigen::VectorXf h;    // diagonal of sqrt(W)*X*(X' V X)^(-1) *X'*sqrt(W)
  Eigen::MatrixXf D;    //  X' V X
  Eigen::MatrixXf covB; // (X' V X)^(-1)
  Eigen::VectorXf beta;
  Eigen::VectorXf delta_beta;

  Eigen::VectorXf y;
  Eigen::VectorXf succ;
  Eigen::VectorXf total;
  // Eigen::VectorXf useless; // keep this, or strange crash bug will happen
};

FirthRegression::FirthRegression()
{
  this->w = new WorkingData;
}

FirthRegression::~FirthRegression()
{
  if (this->w) {
    delete this->w;
    this->w = NULL;
  }
}

Vector & FirthRegression::GetAsyPvalue(){
  int numCov = B.Length();
  pValue.Dimension(B.Length());
  for (int i = 0; i < numCov; i ++){
    double Zstat = B[i] * B[i]/(covB[i][i]);
    // pValue[i] = ndist(Zstat);
    // if (pValue[i] >= 0.5){
    //   pValue[i] = 2*(1-pValue[i]);
    // } else pValue[i] = 2*pValue[i];
    pValue[i] = gsl_cdf_chisq_Q(Zstat, 1.0);
  }
  return(pValue);
}


void FirthRegression::Reset(Matrix& X){
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
  this->w->h.setZero(nr);
  this->w->V.setZero(nr);
  this->w->D.setZero(nc, nc);
  this->w->covB.setZero(nc, nc);
  this->w->beta.setZero(nc);
  this->w->delta_beta.setZero(nc);

  this->w->y.resize(nr);
  this->w->succ.resize(nr);
  this->w->total.resize(nr);
}

bool FirthRegression::FitFirthModel(Matrix & X, Matrix & y, int rnrounds) {
  if (y.cols != 1) {
    fprintf(stderr, "%s:%d Use first column of y\n", __FILE__, __LINE__);
  }
  Vector v(X.rows);
  for (int i = 0; i < X.rows; ++i){
    v[i] = y[i][0];
  }
  return this->FitFirthModel(X, v, rnrounds);
};

bool FirthRegression::FitFirthModel(Matrix & X, Vector & succ, Vector& total, int nrrounds) {
  this-> Reset(X);

  G_to_Eigen(X, &this->w->X);
  G_to_Eigen(succ, &this->w->succ);
  G_to_Eigen(total, &this->w->total);

  int rounds = 0;
  // double lastDeviance, currentDeviance;
  Eigen::MatrixXf xw; // W^(1/2) * X
  // Newton-Raphson
  while (rounds < nrrounds) {
    // beta = beta + solve( t(X)%*%diag(p*(1-p)) %*%X) %*% t(X) %*% (Y-p);
    this->w->eta = this->w->X * this->w->beta;
    this->w->p = (-this->w->eta.array().exp() + 1.0).inverse();
    this->w->V = this->w->p.array() * (1.0 - this->w->p.array()) * this->w->total.array();

    xw = (this->w->V.array().sqrt().matrix().asDiagonal() * this->w->X).eval();
    this->w->D = xw.transpose() * xw; // this->w->X.transpose() * this->w->V.asDiagonal() * this->w->X; // X' V X
    this->w->covB = this->w->D.eval().llt().solve(Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows()));
    // double rel = ((this->w->D * this->w->covB).array() - Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows()).array()).matrix().norm() / this->w->D.rows() / this->w->D.rows();
    // if (rel > 1e-6) { // use relative accuracy to evalute convergence
    if ((this->w->D * this->w->covB - Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows())).norm() > 1e-3) {
      // cannot inverse
      return false;
    }
    this->w->h = (xw.transpose() * this->w->covB * xw).diagonal();
    this->w->r = this->w->X.transpose() * (this->w->succ.array() - this->w->total.array() * this->w->p.array()
                                           + this->w->total.array() * this->w->h.array() * (0.5 - this->w->p.array())).matrix();
    this->w->delta_beta = this->w->covB * this->w->r;
    this->w->beta += this->w->delta_beta;
    if (rounds > 1 && (this->w->beta.norm() > 0 && this->w->delta_beta.norm() / this->w->beta.norm() < 1e-6)) {
      rounds = 0;
      break;
    }
    rounds ++;
  }
  if (rounds == nrrounds) {
    printf("Not enough iterations!");
    return false;
  }

  Eigen_to_G(this->w->beta, &B);
  Eigen_to_G(this->w->covB, &covB);
  Eigen_to_G(this->w->p, &p);
  Eigen_to_G(this->w->V, &V);

  return true;
}

bool FirthRegression::FitFirthModel(Matrix & X, Vector & y, int nrrounds)
{
  this-> Reset(X);

  G_to_Eigen(X, &this->w->X);
  G_to_Eigen(y, &this->w->y);

  int rounds = 0;
  // double lastDeviance, currentDeviance;
  Eigen::MatrixXf xw; // W^(1/2) * X
  // Newton-Raphson  
  while (rounds < nrrounds) {
    // std::cout << "beta = " << this->w->beta << "\n";
    this->w->eta = this->w->X * this->w->beta;
    this->w->p = (1.0 + (-this->w->eta.array()).exp()).inverse();
    this->w->V = this->w->p.array() * (1.0 - this->w->p.array());

    xw = (this->w->V.array().sqrt().matrix().asDiagonal() * this->w->X).eval(); // W^(1/2) * X 
    this->w->D = xw.transpose() * xw; // X' V X
    this->w->covB = this->w->D.eval().ldlt().solve(Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows()));

    // double rel = ((this->w->D * this->w->covB).array() - Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows()).array()).matrix().norm() / this->w->D.rows() / this->w->D.rows();
    // // printf("norm = %g\n", rel);
    // if (rel > 1e-6) { // use relative accuracy to evalute convergence
    if ((this->w->D * this->w->covB - Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows())).norm() > 1e-3) {      
      // cannot inverse
      // printf("Cannot inverse X'WX.\n");
      // printToFile(this->w->D, "matD", "D");
      // printToFile(this->w->covB, "matCovB", "CovB");
      // printToFile(this->w->beta, "matBeta", "beta");      
      return false;
    }
    this->w->h = (xw * this->w->covB * xw.transpose()).diagonal();
    this->w->r = this->w->X.transpose() * (this->w->y -
                                           this->w->p +
                                           (this->w->h.array() * (0.5 - this->w->p.array())).matrix()); // X' (y-mu+h*(1/2-pi))
    this->w->delta_beta = this->w->covB * this->w->r;
    // set dampen coef = 0.5 to reduce the step size
    // otherwise, IWLS may become unstable
    // this seems working well in reality, although more optimal way is to 
    // tune the coef to between 0 and 1.0.
    this->w->beta += this->w->delta_beta * 0.5;
    // printf("norm = %g\n", this->w->delta_beta.norm());
    // use relative accuracy to evalute convergence
    if (rounds > 1 && (this->w->beta.norm() > 0 && this->w->delta_beta.norm() / this->w->beta.norm() < 1e-6)) {
      rounds = 0;
      break;
    }
    rounds ++;
  }
  if (rounds == nrrounds) {
    printf("Not enough iterations!");
    return false;
  }
  // this->w->covB = this->w->D.eval().llt().solve(Eigen::MatrixXf::Identity(this->w->D.rows(), this->w->D.rows()));

  Eigen_to_G(this->w->beta, &B);
  Eigen_to_G(this->w->covB, &covB);
  Eigen_to_G(this->w->p, &p);
  Eigen_to_G(this->w->V, &V);

  return true;
}



#if 0
double FirthRegression::GetDeviance() {
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

double FirthRegression::GetDeviance(Matrix & X, Vector & y)
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

double FirthRegression::GetDeviance(Matrix & X, Vector & succ, Vector& total)
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
#endif
