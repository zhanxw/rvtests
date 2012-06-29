#include "Skat.h"

#include "MathMatrix.h"
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

#define ZBOUND 1e-30

void Eigen_to_G(Eigen::MatrixXf &EigenM, Matrix* GM);
void Eigen_to_G(Eigen::VectorXf &EigenV, Vector* GV);
void G_to_Eigen(Vector &GV, Eigen::VectorXf* EigenV);
void G_to_Eigen(Matrix &GM, Eigen::MatrixXf* EigenM);

void MatrixSqrt(Eigen::MatrixXf& in, Eigen::MatrixXf* out);


int Skat::Fit(Vector & res_G,   // residual under NULL -- may change when permuting
              Vector& v_G,      // variance under NULL -- may change when permuting
              Matrix& X_G,      // covariance
              Matrix & G_G,     // genotype
              Vector &w_G)     // weight 
{
  this->nPeople = X_G.rows;
  this->nMarker = G_G.cols;
  this->nCovariate = X_G.cols;

  // calculation w_sqrt
  G_to_Eigen(w_G, &this->w_sqrt);
  w_sqrt.cwiseSqrt();

  // calculate K
  Eigen::MatrixXf G;
  G_to_Eigen(G_G, &G);
  this->K_sqrt.noalias() = w_sqrt.asDiagonal() * G.transpose();

  // calculate Q
  Eigen::VectorXf res;
  G_to_Eigen(res_G, &res);
  this->Q = (this->K_sqrt * res).squaredNorm();
  
  // calculate P0
  Eigen::VectorXf v;
  G_to_Eigen(v_G, &v);
  if (this->nCovariate == 1) {
    P0 = v.asDiagonal();
    P0 -= v * v.transpose() / v.sum();
  } else {
    Eigen::MatrixXf X;
    G_to_Eigen(X_G, &X);    
    Eigen::MatrixXf XtV ;        // X^t V
    XtV.noalias() = X.transpose() * v.asDiagonal();
    P0 = v.asDiagonal();
    P0 -= XtV.transpose() * ( ( XtV * X ).inverse() ) * XtV;
  }
  // eigen decomposition
  es.compute( K_sqrt * P0 * K_sqrt.transpose());
  this->mixChiSq.reset();
  int r_ub = std::min(nPeople, nMarker);
  int r = 0; // es.eigenvalues().size();
  int eigen_len = es.eigenvalues().size();
  for(int i=eigen_len-1; i>=0; i--)
  {
    if (es.eigenvalues()[i] > ZBOUND && r<r_ub) {
      mixChiSq.addLambda(es.eigenvalues()[i]);
      r++;
    }
	else break;
  }
  // calculate p-value
  this->pValue = this->mixChiSq.getPvalue(this->Q);
  return 0;
};

double Skat::GetQFromNewResidual(Vector & res_G){   // residual under NULL -- may change when permuting
  // calculate Q
  Eigen::VectorXf res;
  G_to_Eigen(res_G, &res);
  this->Q = (this->K_sqrt * res).squaredNorm();

  // calculate p-value
  this->pValue = this->mixChiSq.getPvalue(this->Q);
};

#if 0
int Skat::CalculatePValue(Vector & y_G, Vector& y0_G, Matrix& X_G, Vector& v_G,
                          Matrix & G_G, Vector &w_G) {
  // Eigen::VectorXf y;
  // Eigen::VectorXf y0;
  Eigen::MatrixXf X;
  Eigen::VectorXf v;
  Eigen::MatrixXf G;
  Eigen::VectorXf w;

  // G_to_Eigen(y_G, y);
  // G_to_Eigen(y0_G, y0);
  G_to_Eigen(X_G, X);
  G_to_Eigen(v_G, v);
  G_to_Eigen(G_G, G);
  G_to_Eigen(w_G, w);

  // get residual
  int yLen = y_G.Length();
  res.resize(yLen);
  for (int i = 0; i < yLen; i++) {
    res(i) = y_G[i] - y0_G[i];
  }

  // cache P0
  // P0 = V - V * X * (X' * V * X)^{-1} * X' * V        NOTE: V is symmetric
  // when X is column vector of 1,
  // P0 = V - V * X * X' V / \sum(v_i)        v_i is ith element on V
  //    = V - v * v' / \sum(v_i)              v is column vector of v_i
  if (X_G.cols == 1) {
    P0 = v.asDiagonal();
    P0 -= v * v.transpose() / v.sum();
  } else {
    Eigen::MatrixXf XtV;        // X^t V
    XtV = X.transpose() * v.asDiagonal();
    P0 = v.asDiagonal();
    P0 -= XtV.transpose() * ( ( XtV * X ).inverse() ) * XtV;
  }
  //outputMat(P0, "mat.P0");

  // prep parameters to qf()
  lambda = new double[v.size()];
  noncen = new double[v.size()];
  df = new int[v.size()];
  lambda_size = v.size();

  //Re-allocate memory if # of markers is > # of samples
  if(lambda_size < G.cols())
  {
	lambda_size = G.cols();
	delete [] lambda; delete [] noncen; delete [] df;
    lambda = new double[lambda_size];
    noncen = new double[lambda_size];
    df = new int[lambda_size];
  }
  // w_sqrt <- W
  w_sqrt.resize(w_G.Length());
  for (int i = 0; i < w_G.Length(); i++) {
    if (w_G[i] > ZBOUND)
      w_sqrt(i) = sqrt(w_G[i]);
    else
      w_sqrt(i) = 0.0;
  }

  // get K_sqrt
  K_sqrt = w_sqrt.asDiagonal() * G.transpose();

  // get Q = (K_sqrt * res)' * (K_sqrt * res)
  Eigen::VectorXf tmp = (K_sqrt * res);
  this->Q = tmp.squaredNorm();

  // get P1 = w
  Eigen::MatrixXf P1 = w_sqrt.asDiagonal() * G.transpose() * P0 * G * w_sqrt.asDiagonal();

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;
  es.compute(P1);

  // Input to qf
  int r_ub = G.rows() > G.cols() ? G.cols() : G.rows();
  int r = 0; // es.eigenvalues().size();
  int eigen_len = es.eigenvalues().size();
  for(int i=eigen_len-1; i>=0; i--)
  {
    if (es.eigenvalues()[i] > ZBOUND && r<r_ub) {
      lambda[r] = es.eigenvalues()[i];
      noncen[r] = 0.0;
      df[r] = 1;
      r++;
    }
	else break;
  }


  // Output from qf
  double pvalue;
  int fault;
  double trace[7];

  // note: qf give distribution but we want p-value.
  this->pValue = 1.0 - qf(lambda, noncen, df, r, sigma, Q, lim, acc, trace, &fault);
  if(this->pValue>1.0) this->pValue = 1.0; //this occurs when eigen values are very large

  // TIME();

  if (fault) {
    return fault;
  }

  return 0;
}
#endif

void G_to_Eigen(Matrix& GM, Eigen::MatrixXf* _EigenM)
{
  Eigen::MatrixXf& EigenM = *_EigenM;
  EigenM.resize(GM.rows, GM.cols);
  for(int i=0; i<GM.rows; i++)
    for(int j=0; j<GM.cols; j++)
      EigenM(i, j) = GM[i][j];
}
void G_to_Eigen(Vector& GV, Eigen::VectorXf* _EigenV)
{
  Eigen::VectorXf& EigenV = *_EigenV;
  EigenV.resize(GV.Length());
  for(int i=0; i<GV.Length(); i++)
    EigenV(i) = GV[i];
}

void Eigen_to_G(Eigen::MatrixXf& EigenM, Matrix* _GM)
{
  Matrix& GM = *_GM;
  GM.Dimension(EigenM.rows(), EigenM.cols());
  for(int i=0; i<GM.rows; i++)
    for(int j=0; j<GM.cols; j++)
      GM[i][j] = EigenM(i, j);
}


void Eigen_to_G(Eigen::VectorXf& EigenV, Vector* _GV)
{
  Vector& GV = *_GV;
  EigenV.resize(GV.Length());
  for(int i=0; i<GV.Length(); i++)
    EigenV(i) = GV[i];
}

/**
 *     // NOTE: since @param may not be full rank, we will use SVD
 */
void MatrixSqrt(Eigen::MatrixXf& in, Eigen::MatrixXf* _out) {
  Eigen::MatrixXf& out = *_out;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(in);
  out = es.operatorSqrt();
  // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(in);
  // Eigen::VectorXf d = es.eigenvalues();;
  // for (int i = 0; i < d.size(); i++){
  //   if (d(i) > ZBOUND) {
  //     d(i) = sqrt(d[i]);
  //   }else{
  //     d(i) = 0.0;
  //   }
  // }
  // out = es.eigenvectors() * d.asDiagonal() *  es.eigenvectors().transpose();
}
