#include "FamSkat.h"
#include "MathMatrix.h"
#include "EigenMatrix.h"
#include "EigenMatrixInterface.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
//#include <Eigen/Dense>
#include "MixtureChiSquare.h"
#include "FastLMM.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#undef DEBUG
// #define DEBUG
#ifdef DEBUG
#include <fstream>
#endif

#define ZBOUND 1e-30

class FamSkat::FamSkatImpl{
 public:
  FamSkatImpl() :
      lmm(FastLMM::SCORE, FastLMM::MLE),
      beta1(1), beta2(25) {
  }
  int FitNullModel(Matrix& Xnull, Matrix& y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
    return this->lmm.FitNullModel(Xnull, y, kinshipU, kinshipS);
  }
  int TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol, Vector& weight,
                    const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
    this->nPeople = Xnull.rows;
    this->nMarker = Xcol.cols;
    this->nCovariate = Xnull.cols;

    // set up weight
    G_to_Eigen(weight, &this->weight);
    setupWeight(kinshipU, kinshipS, Xcol);

    G_to_Eigen(Xnull, &this->X);
    G_to_Eigen(Xcol,  &this->G);
    G_to_Eigen(y,     &this->yMat);

    const double sigma2 = lmm.GetSigmaG2();
    EigenMatrix tmp;
    lmm.GetBeta(&tmp);
    this->beta = tmp.mat;
    const double delta = lmm.GetDelta();
    const Eigen::MatrixXf& U =  kinshipU.mat;
    const Eigen::MatrixXf& lambda =  kinshipS.mat;
    // Sigma = U * (S + delta) * U' * sigma2
    // Q = || w^{1/2} * G' * Sigma^{-1} * (y - X * beta) ||
    //   = || w^{1/2} * G' * U * (S + delta)^{-1} * U' * (y - X * beta) /sigma2 ||
    Q = (this->weight.asDiagonal() *
         G.transpose() *
         U *
         (lambda.array() + delta).inverse().matrix().asDiagonal()  *
         U.transpose() *
         (yMat - X * beta) /
         sigma2
         ).col(0).squaredNorm();

    // P0 = Sigma - X * (X' * Sigma^{-1} * X)^{-1} * X'
    Eigen::MatrixXf Sigma = U *
                            (lambda.array() + delta).matrix().asDiagonal() *
                            U.transpose() * sigma2;
    Eigen::MatrixXf SigmaInv = U *
                               (lambda.array() + delta).inverse().matrix().asDiagonal() *
                               U.transpose() / sigma2;

    Eigen::MatrixXf k = this->weight.asDiagonal() *
                        G.transpose() *
                        SigmaInv;
    Eigen::MatrixXf P0 = Sigma - X * (X.transpose() * SigmaInv * X).inverse() * X.transpose();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;
    es.compute( k * P0 * k.transpose());

#ifdef DEBUG
    std::ofstream kin("kin");
    kin << U * lambda.asDiagonal() * U.transpose();
    kin << "\n";
    kin <<       U *
        (lambda.array() + delta).inverse().matrix().asDiagonal()  *
        U.transpose();

    kin.close();


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

  }
#if 0
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
#endif

  // double GetQFromNewResidual(Vector & res_G)   // e.g. permuted residual under NULL
  // {
  //   // calculate Q
  //   Eigen::VectorXf res;
  //   G_to_Eigen(res_G, &res);
  //   double Q = (this->K_sqrt * res).squaredNorm();

  //   return Q;
  // }

  double GetPvalue() const {return this->pValue;};

  double GetQ() const {return this->Q;};
 private:
  void setupWeight(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                   Matrix& geno) {
    const int p = geno.cols;
    this->weight.resize(p);
    for (int i = 0; i < p; ++i) {
      double af = lmm.FastGetAF(kinshipU, kinshipS, geno, i);
      this->weight[i] =  gsl_ran_beta_pdf(af, this->beta1, this->beta2);  /// default SKAT use beta(MAF, 1, 25)
    }
  }
 private:
  Eigen::VectorXf weight;     // p by 1 matrix (p: number of variant) NOTE: not squared as the paper notation
  Eigen::MatrixXf X;
  Eigen::MatrixXf G;
  Eigen::MatrixXf yMat;
  Eigen::MatrixXf beta;

  //Eigen::MatrixXf K;        // G * W * G'
  // Eigen::MatrixXf K_sqrt;     // W^{0.5} * G' ----> K = K_sqrt' * K_sqrt
  Eigen::VectorXf w_sqrt;     // W^{0.5} * G' ----> K = K_sqrt' * K_sqrt
  Eigen::MatrixXf P0;         // V - VX ( X' V X)^{-1} X V
  Eigen::VectorXf res;        // residual

  FastLMM lmm;

  double beta1;
  double beta2;

  int nPeople;
  int nMarker;
  int nCovariate;
  MixtureChiSquare mixChiSq;

  double pValue;
  double Q;

};

//////////////////////////////////////////////////
// FamSkat class
FamSkat::FamSkat() {
  this->skatImpl =  new FamSkatImpl;
}
FamSkat::~FamSkat() {
  delete this->skatImpl;
}
int FamSkat::FitNullModel(Matrix& Xnull, Matrix& y,
                          const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
  return this->skatImpl->FitNullModel(Xnull, y, kinshipU, kinshipS);
}
int FamSkat::TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol, Vector& weight,
                           const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
  return this->skatImpl->TestCovariate(Xnull, y, Xcol, weight, kinshipU, kinshipS);
}

// void FamSkat::Reset() {
//   this->skatImpl->Reset();
// }

// double FamSkat::GetQFromPermutation();
//   return this->skatImpl->GetQFromNewResidual(res_G);
// }

double FamSkat::GetPvalue() const //  {return this->pValue;};
{
  return this->skatImpl->GetPvalue();
}

double FamSkat::GetQ() const // {return this->Q;};
{
  return this->skatImpl->GetQ();
}
