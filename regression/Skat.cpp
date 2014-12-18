#include "Skat.h"
#include "MathMatrix.h"
#include "EigenMatrixInterface.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
//#include <Eigen/Dense>
#include "MixtureChiSquare.h"

// #define DEBUG
#undef DEBUG
#ifdef DEBUG
#include <fstream>
#endif

#define ZBOUND 1e-30
void MatrixSqrt(Eigen::MatrixXf& in, Eigen::MatrixXf* out);

class Skat::SkatImpl{
 public:
  void Reset() {
    this->pValue = -999.0;
  };
  int Fit(Vector & res_G,   // residual under NULL -- may change when permuting
          Vector& v_G,      // variance under NULL -- may change when permuting
          Matrix& X_G,      // covariance
          Matrix & G_G,     // genotype
          Vector &w_G)      // weight
  {
    this->nPeople = X_G.rows;
    this->nMarker = G_G.cols;
    this->nCovariate = X_G.cols;

    // calculation w_sqrt
    G_to_Eigen(w_G, &this->w_sqrt);
    w_sqrt = w_sqrt.cwiseSqrt();

    // calculate K = G * W * G'
    Eigen::MatrixXf G;
    G_to_Eigen(G_G, &G);
    this->K_sqrt.noalias() = w_sqrt.asDiagonal() * G.transpose();

    // calculate Q = ||res * K||
    Eigen::VectorXf res;
    G_to_Eigen(res_G, &res);
    this->Q = (this->K_sqrt * res).squaredNorm();

    // calculate P0 = V - V X (X' V X)^(-1) X' V
    Eigen::VectorXf v;
    G_to_Eigen(v_G, &v);
    if (this->nCovariate == 1) {
      P0 = - v * v.transpose() / v.sum();
      // printf("dim(P0) = %d, %d\n", P0.rows(), P0.cols());
      // printf("dim(v) = %d\n", v.size());
      P0.diagonal() += v;
      // printf("dim(v) = %d\n", v.size());
    } else {
      Eigen::MatrixXf X;
      G_to_Eigen(X_G, &X);
      Eigen::MatrixXf XtV ;        // X^t V
      XtV.noalias() = X.transpose() * v.asDiagonal();
      P0 = - XtV.transpose() * ( ( XtV * X ).inverse() ) * XtV;
      P0.diagonal() += v;
    }
    // dump();
    // Eigen::MatrixXf tmp = K_sqrt * P0 * K_sqrt.transpose();
    // dumpToFile(tmp, "out.tmp");
    // eigen decomposition
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;
    es.compute( K_sqrt * P0 * K_sqrt.transpose());

#ifdef DEBUG
    std::ofstream k("K");
    k << K_sqrt;
    k.close();
#endif
    // std::ofstream p("P0");
    // p << P0;
    // p.close();

    this->mixChiSq.reset();
    int r_ub = std::min(nPeople, nMarker);
    int r = 0; // es.eigenvalues().size();
    int eigen_len = es.eigenvalues().size();
    for(int i=eigen_len-1; i>=0; i--)
    {
      if (es.eigenvalues()[i] > ZBOUND && r<r_ub) {
        this->mixChiSq.addLambda(es.eigenvalues()[i]);
        r++;
      }
      else break;
    }
    // calculate p-value
    this->pValue = this->mixChiSq.getPvalue(this->Q);
    return 0;
  };


  double GetQFromNewResidual(Vector & res_G)   // e.g. permuted residual under NULL
  {
    // calculate Q
    Eigen::VectorXf res;
    G_to_Eigen(res_G, &res);
    double Q = (this->K_sqrt * res).squaredNorm();

    return Q;
  }

  double GetPvalue() const {return this->pValue;};

  double GetQ() const {return this->Q;};

 private:
  //Eigen::MatrixXf K;        // G * W * G'
  Eigen::MatrixXf K_sqrt;     // W^{0.5} * G' ----> K = K_sqrt' * K_sqrt
  Eigen::VectorXf w_sqrt;     // W^{0.5} 
  Eigen::MatrixXf P0;         // V - VX ( X' V X)^{-1} X V
  Eigen::VectorXf res;        // residual

  int nPeople;
  int nMarker;
  int nCovariate;
  MixtureChiSquare mixChiSq;

  double pValue;
  double Q;
};
Skat::Skat() {
  this->skatImpl =  new SkatImpl;
}
Skat::~Skat() {
  delete this->skatImpl;
}
void Skat::Reset() {
  this->skatImpl->Reset();
}

int Skat::Fit(Vector & res_G,   // residual under NULL -- may change when permuting
              Vector& v_G,      // variance under NULL -- may change when permuting
              Matrix& X_G,      // covariance
              Matrix & G_G,     // genotype
              Vector &w_G)     // weight
{
  return this->skatImpl->Fit(res_G, v_G, X_G, G_G, w_G);
};

double Skat::GetQFromNewResidual(Vector & res_G){   // residual under NULL -- may change when permuting
  return this->skatImpl->GetQFromNewResidual(res_G);
}
double Skat::GetPvalue() const //  {return this->pValue;};
{
  return this->skatImpl->GetPvalue();
}

double Skat::GetQ() const // {return this->Q;};
{
  return this->skatImpl->GetQ();
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
