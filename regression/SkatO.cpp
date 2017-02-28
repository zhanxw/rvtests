#include "SkatO.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <vector>

#include "gsl/gsl_cdf.h"      // use gsl_cdf_chisq_Q
#include "gsl/gsl_randist.h"  // gsl_ran_chisq_pdf

#include "EigenMatrixInterface.h"
#include "GSLIntegration.h"
#include "MathMatrix.h"
#include "MixtureChiSquare.h"

// #define DEBUG
#undef DEBUG
#ifdef DEBUG
#include <fstream>
#include "MatrixOperation.h"
template <class T>
void dumpEigen(const char* fn, T& m) {
  std::ofstream k(fn);
  k << m;
  k.close();
}
template <class T>
void printEigen(T& m) {
  std::ofstream k;
  k << m;
}
#endif

#define ZBOUND 1e-30

struct Moment {
  double muQ;
  double varQ;
  double df;
};

double integrandDavies(double x, void* param);
double integrandLiu(double x, void* param);

class SkatO::SkatOImpl {
 public:
  SkatOImpl() {
    for (int i = 0; i <= 10; ++i) {
      rhosOriginal.push_back(1.0 * i / 10);
    }
    this->nRho = rhosOriginal.size();
  }

  void Reset() { this->pValue = -999.0; };
  void setBinaryOutcome() { this->binaryOutcome = true; }
  void setQuantitativeOutcome() { this->binaryOutcome = false; }
  bool isBinary() { return binaryOutcome; }

  int FitSKAT(Eigen::VectorXd res,  // residual under NULL
              Eigen::VectorXd v,    // variance under NULL
              Eigen::MatrixXd X,    // covariance
              Eigen::MatrixXd G,    // genotype
              Eigen::VectorXd w)    // weight
  {
    this->rho = 0;

    G = G * w.asDiagonal();
    Eigen::MatrixXd temp = res.transpose() * G;
    this->Q = (temp * temp.transpose())(0, 0);

    if (!isBinary()) {
      double s2 = res.squaredNorm() / (nPeople - 1);
      this->Q /= s2;
    }
    this->Q /= 2.;

    // SKAT R:
    // W.1 = t(Z) %*% Z - (t(Z) %*% X1)%*%solve(t(X1) %*% X1) %*% *(t(X1) %*% Z)
    // # t (Z) P0 Z
    // calculate Z1
    if (!isBinary()) {
      W = G.transpose() * G -
          (G.transpose() * X) *
              (X.transpose() * X).ldlt().solve(X.transpose() * G);
    } else {
      W = G.transpose() * v.asDiagonal() * G -
          (G.transpose() * v.asDiagonal() * X) *
              (X.transpose() * v.asDiagonal() * X)
                  .ldlt()
                  .solve(X.transpose() * v.asDiagonal() * G);
    }
    W = W / 2;  // follow SKAT R package convension to divide 2 here.
    if (getEigen(W, &lambda)) {  // error can occur when lambda are all zeros
      return -1;
    }
    this->pValue = getPvalDavies(Q, lambda);
    return 0;
  }

  int Fit(Vector& res_G,  // residual under NULL
          Vector& v_G,    // variance under NULL
          Matrix& X_G,    // covariance
          Matrix& G_G,    // genotype
          Vector& w_G)    // weight
  {
    this->nPeople = X_G.rows;
    this->nMarker = G_G.cols;
    this->nCovariate = X_G.cols;

    // calculation w_sqrt
    G_to_Eigen(res_G, &this->res);
    G_to_Eigen(v_G, &this->v);
    G_to_Eigen(X_G, &this->X);
    G_to_Eigen(G_G, &this->G);
    G_to_Eigen(w_G, &this->w);

    if (G.cols() == 1) {
      return FitSKAT(res, v, X, G, w);
    }
    // avoid rho = 1.0
    capRhos();

    G = G * w.asDiagonal();

    // set up R_rho
    R_rhos.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      getRrho(rhos[i], &R_rhos[i], G.cols());
    }

    double s2;
    if (isBinary()) {
      s2 = 1;
    } else {  // continuous trait
      s2 = res.norm();
      s2 = (s2 * s2) / (nPeople - 1);
    }

    // calculate Q
    Qs.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      Eigen::MatrixXd v = res.transpose() * G;
      Qs[i] = (v * R_rhos[i] * v.transpose())(0, 0);
      Qs[i] /= s2;
      Qs[i] /= 2.0;  // follow SKAT R package convention (divides 2)
    }

    // calculate Z1
    if (!isBinary()) {
      Z1 = G - X * (X.transpose() * X).ldlt().solve(X.transpose() * G);
    } else {
      Eigen::VectorXd v_sqrt = v.cwiseSqrt();
      Z1 = v_sqrt.asDiagonal() * G -
           v_sqrt.asDiagonal() * X *
               (X.transpose() * v.asDiagonal() * X)
                   .ldlt()
                   .solve(X.transpose() * v.asDiagonal() * G);
    }
    Z1 = Z1 / sqrt(2);  // follow SKAT R package convention (divides sqrt{2})

    // calculate labmda.rho
    lambdas.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      Eigen::LLT<Eigen::MatrixXd> chol;
      chol.compute(R_rhos[i]);  // R_rhos[i] = L * L^T
      Eigen::MatrixXd Z2 = Z1 * chol.matrixL();
      Eigen::MatrixXd K = Z2.transpose() * Z2;
      if (getEigen(K, &lambdas[i])) {
        // error occured,
        // e.g. G is in the column space of Z => Z1 = 0 => K is all zeros
        //      this can happen when many covariates are used
        return -1;
      }
    }

    // calculate some parameters (for Z(I-M)Z part)
    Eigen::MatrixXd z_bar = Z1.rowwise().sum() / Z1.cols();
    double z_norm = z_bar.squaredNorm();
    Eigen::MatrixXd ZMZ = z_bar.transpose() * Z1;
    ZMZ = (ZMZ.transpose() * ZMZ) / z_norm;
    Eigen::MatrixXd ZIMZ = Z1.transpose() * Z1 - ZMZ;  // Z(I-M)Z' = ZZ' - ZMZ'
    if (getEigen(ZIMZ, &this->lambda)) {
      return -1;
    }
    this->VarZeta = 4.0 * (ZMZ.array() * ZIMZ.array()).sum();
    this->MuQ = lambda.sum();
    this->VarQ = 2.0 * (lambda.array() * lambda.array()).sum() + VarZeta;

    double temp = (lambda.array() * lambda.array()).sum();
    double KerQ =
        (lambda.array() * lambda.array() * lambda.array() * lambda.array())
            .sum() /
        temp / temp * 12;
    this->Df = 12 / KerQ;

    // calculate tau
    taus.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      taus[i] = nMarker * nMarker * rhos[i] * z_norm +
                (1.0 - rhos[i]) *
                    (z_bar.transpose() * Z1).array().square().sum() / z_norm;
    }

    // calculate moments
    moments.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      getMoment(lambdas[i], &moments[i]);
    }

    // calculate p for each rho
    pvals.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      pvals[i] = getPvalByMoment(Qs[i], moments[i]);
    }

    // calculat min p
    double minP = pvals[0];
    int minIndex = 0;
    for (int i = 1; i < nRho; ++i) {
      if (pvals[i] < minP) {
        minP = pvals[i];
        minIndex = i;
      }
    }
    this->rho = rhos[minIndex];
    this->Q = Qs[minIndex];

    // calculate Q_minP
    Qs_minP.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      Qs_minP[i] = getQvalByMoment(minP, moments[i]);
    }

    // integrate
    Integration integration;
    integration.setEpsAbs(1e-25);
    integration.setEpsRel(
        0.0001220703);  // this is the default value of epsrel in R
                        // .Machine$double.eps^0.25
    gsl_function F;
    F.function = integrandDavies;
    F.params = this;
    if (integration.integrateLU(F, 0., 40.)) {
#ifdef DEBUG
      fprintf(stderr, "%s:%d integration failed\n", __FILE__, __LINE__);
#endif
      F.function = integrandLiu;
      F.params = this;
      if (integration.integrateLU(F, 0., 40.)) {
#ifdef DEBUG
        fprintf(stderr, "%s:%d integration failed\n", __FILE__, __LINE__);
#endif
      }
    };
    this->pValue = 1.0 - integration.getResult();

    // verify p-values
    // SKAT R: "Since SKAT-O is between burden and SKAT, SKAT-O p-value should
    // be <= min(p-values) * 2
    // To correct conservatively, we use min(p-values) * 3"
    const int multi = (nRho < 3) ? 2 : 3;
    if (nRho)
      if (this->pValue <= 0) {
        double p = minP * multi;
        if (pValue < p) {
          pValue = p;
        }
      }
    if (this->pValue == 0.0) {  // for some i, pvals[i] == 0.0
      pValue = pvals[0];
      for (int i = 1; i < nRho; ++i) {
        if (pvals[i] > 0 && pvals[i] < pValue) {
          pValue = pvals[i];
        }
      }
    }

    uncapRhos();
    return 0;
  }

  double getPvalDavies(double Q, const Eigen::MatrixXd& lambda) {
    double pval = 0.0;
    this->mixChiSq.reset();
    for (int i = 0; i < lambda.rows(); ++i) {
      this->mixChiSq.addLambda(lambda(i, 0));
    }
    pval = this->mixChiSq.getPvalue(Q);
    return pval;
  }

  double getPvalLiu(double Q, const Eigen::MatrixXd& lambda) {
    double pval = 0.0;
    this->mixChiSq.reset();
    for (int i = 0; i < lambda.rows(); ++i) {
      this->mixChiSq.addLambda(lambda(i, 0));
    }
    pval = this->mixChiSq.getLiuPvalue(Q);
    return pval;
  }

  double computeIntegrandDavies(double x) {
    double kappa = DBL_MAX;
    for (int i = 0; i < nRho; ++i) {
      double v = (Qs_minP[i] - taus[i] * x) / (1.0 - rhos[i]);
      if (i == 0) {
        kappa = v;
      }
      if (v < kappa) {
        kappa = v;
      }
    }
    double temp;
    if (kappa > lambda.sum() * 10000) {
      temp = 0.0;
    } else {
      double Q = (kappa - MuQ) * sqrt(VarQ - VarZeta) / sqrt(VarQ) + MuQ;
      temp = getPvalDavies(Q, lambda);
      if (temp <= 0.0 || temp == 1.0) {
        temp = getPvalLiu(Q, lambda);
      }
    }
    return (1.0 - temp) * gsl_ran_chisq_pdf(x, 1.0);
  }

  double computeIntegrandLiu(double x) {
    double kappa = DBL_MAX;
    for (int i = 0; i < nRho; ++i) {
      double v = (Qs_minP[i] - taus[i] * x) / (1.0 - rhos[i]);
      if (v < kappa) {
        kappa = v;
      }
    }
    double Q = (kappa - MuQ) / sqrt(VarQ) * sqrt(2.0 * Df) + Df;
    return (gsl_cdf_chisq_P(Q, Df) * gsl_ran_chisq_pdf(x, 1.0));
  }
  double GetPvalue() const { return this->pValue; };
  double GetQ() const { return this->Q; };
  double GetRho() const { return this->rho; }

 private:
  int getRrho(double rho, Eigen::MatrixXd* R, int dim) {
    (*R).setConstant(dim, dim, rho);
    for (int i = 0; i < dim; ++i) {
      (*R)(i, i) = 1.0;
    }
    return 0;
  }
  int getEigen(Eigen::MatrixXd& k, Eigen::MatrixXd* lambda) {
    es.compute(k);
    Eigen::VectorXd values = es.eigenvalues();
    int n = values.size();
    int numNonZero = 0;
    double sumNonZero = 0.;
    for (int i = 0; i < n; ++i) {
      if (values(i) > 0) {
        ++numNonZero;
        sumNonZero += values(i);
      }
    }
    if (numNonZero == 0) {
      return -1;
    }
    double t = sumNonZero / numNonZero / 100000;
    int numKeep = n;
    for (int i = 0; i < n; ++i) {
      // values are in increase order, so we can stop early
      if (values[i] < t) {
        --numKeep;
      } else {
        break;
      }
    }

    (*lambda).resize(numKeep, 1);
    for (int i = 0; i < numKeep; ++i) {
      (*lambda)(i, 0) = values(n - 1 - i);
    }

    return 0;
  }
  int getMoment(const Eigen::MatrixXd& lambda, Moment* moment) {
    Moment& m = *moment;
    double c[4];
    Eigen::ArrayXd la = lambda.array();
    c[0] = la.sum();
    c[1] = (la.square()).sum();
    c[2] = (la.square() * la).sum();
    c[3] = (la.square().square()).sum();

    m.muQ = c[0];
    double sigmaQ = sqrt(2 * c[1]);
    double s1 = c[2] / c[1] / sqrt(c[1]);
    double s2 = c[3] / (c[1] * c[1]);

    // double beta1 = sqrt(8.)*s1;
    // double beta2 = 12*s2;

    double a, d, l;
    if (s1 * s1 > s2) {
      a = 1 / (s1 - sqrt(s1 * s1 - s2));
      d = (s1 * a - 1.0 * a * a);
      l = a * a - 2 * d;
    } else {
      // type1 = 1;
      l = 1. / s2;
      a = sqrt(l);
      d = 0;
    }
    // double muX =l+d;
    // double sigmaX=sqrt(2.) *a;
    m.varQ = sigmaQ * sigmaQ;
    m.df = l;
    return 0;
  }

  double getPvalByMoment(double Q, const Moment& moment) {
    double muQ = moment.muQ;
    double varQ = moment.varQ;
    double df = moment.df;
    double Q_Norm = (Q - muQ) / sqrt(varQ) * sqrt(2. * df) + df;
    double pval = gsl_cdf_chisq_Q(Q_Norm, df);
    return pval;
  }

  double getQvalByMoment(double min_pval, const Moment& moment) {
    double muQ = moment.muQ;
    double varQ = moment.varQ;
    double df = moment.df;
    double q_org = gsl_cdf_chisq_Qinv(min_pval, df);
    double q = (q_org - df) / sqrt(2. * df) * sqrt(varQ) + muQ;
    return q;
  }

  void capRhos() {
    // cap max(rho) to 0.999 to avoid rank-deficiency
    rhos.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      if (rhosOriginal[i] > 0.999) {
        rhos[i] = 0.999;
      } else {
        rhos[i] = rhosOriginal[i];
      }
    }
  }
  void uncapRhos() {
    // cap max(rho) to 0.999 to avoid rank-deficiency
    for (int i = 0; i < nRho; ++i) {
      rhos[i] = rhosOriginal[i];
    }
    if (this->rho >= 0.999) {
      this->rho = 1.;
    }
  }

 private:
  bool binaryOutcome;

  Eigen::VectorXd res;  // residual
  Eigen::VectorXd v;
  Eigen::MatrixXd X;
  Eigen::MatrixXd G;
  Eigen::VectorXd w;  // residual

  Eigen::MatrixXd W;  // for SKAT calculation
  Eigen::MatrixXd Z1;
  std::vector<double> rhos;
  std::vector<double> rhosOriginal;
  std::vector<Eigen::MatrixXd> R_rhos;
  std::vector<Eigen::MatrixXd> lambdas;
  std::vector<Moment> moments;
  std::vector<double> Qs;
  std::vector<double> Qs_minP;
  std::vector<double> taus;
  std::vector<double> pvals;

  int nPeople;
  int nMarker;
  int nCovariate;
  int nRho;

  MixtureChiSquare mixChiSq;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;

  double MuQ;
  double VarQ;
  double VarZeta;
  double Df;
  Eigen::MatrixXd lambda;

  double pValue;
  double Q;
  double rho;
};
SkatO::SkatO() { this->skatoImpl = new SkatOImpl; }
SkatO::~SkatO() { delete this->skatoImpl; }
void SkatO::Reset() { this->skatoImpl->Reset(); }

int SkatO::Fit(Vector& res_G,    // residual under NULL
               Vector& v_G,      // variance under NULL
               Matrix& X_G,      // covariance
               Matrix& G_G,      // genotype
               Vector& w_G,      // weight
               const char* type  // response type
               ) {
  if (!type) return -1;
  switch (type[0]) {
    case 'D':
      this->skatoImpl->setBinaryOutcome();
      break;
    case 'C':
      this->skatoImpl->setQuantitativeOutcome();
      break;
    default:
      return -1;
  }
  return this->skatoImpl->Fit(res_G, v_G, X_G, G_G, w_G);
};

double SkatO::GetPvalue() const { return this->skatoImpl->GetPvalue(); }

double SkatO::GetQ() const { return this->skatoImpl->GetQ(); }

double SkatO::GetRho() const { return this->skatoImpl->GetRho(); }

double integrandDavies(double x, void* param) {
  SkatO::SkatOImpl* p = (SkatO::SkatOImpl*)param;
  return p->computeIntegrandDavies(x);
}
double integrandLiu(double x, void* param) {
  SkatO::SkatOImpl* p = (SkatO::SkatOImpl*)param;
  return p->computeIntegrandLiu(x);
}
