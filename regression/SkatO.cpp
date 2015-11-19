#include "SkatO.h"

#include <Eigen/Core>
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>
#include <vector>

#include "gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q
#include "gsl/gsl_randist.h" // gsl_ran_chisq_pdf

#include "EigenMatrixInterface.h"
#include "GSLIntegration.h"
#include "MathMatrix.h"
#include "MixtureChiSquare.h"


#define DEBUG
// #undef DEBUG
#ifdef DEBUG
#include <fstream>
template <class T>
void dumpEigen(const char* fn, T& m) {
  std::ofstream k(fn);
  k << m;
  k.close();
}
template <class T>
void printEigen(T& m) {
  std::ostream k;
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
    for(int i = 0; i <= 10; ++i) {
      rhos.push_back(1.0 * i / 10);
    }
    this->nRho = rhos.size();
    // cap max(rho) to 0.999 to avoid rank-deficiency
    for(int i = 0; i < nRho; ++i) {
      if (rhos[i] > 0.999) {
        rhos[i] = 0.999;
      }
    }
  }
  void Reset() { this->pValue = -999.0; };
  int Fit(Vector& res_G,  // residual under NULL -- may change when permuting
          Vector& v_G,    // variance under NULL -- may change when permuting
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

    G = G * w.asDiagonal();

    // set up R_rho
    R_rhos.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      getRrho(rhos[i], &R_rhos[i], G.cols());
    }

    float s2;
    if (isBinary()) {
      s2 = 1;
    } else { // continuous trait
      s2 = res.norm();
      s2 = (s2 * s2) / (nPeople - 1);
    }

    // calculate Q
    Qs.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      Eigen::MatrixXf v = res.transpose() * G;
      Qs[i] = (v * R_rhos[i] * v.transpose())(0, 0);
      Qs[i] /= s2;
      Qs[i] /= 2.0; // follow SKAT R package convention (divides 2)
    }

    // calculate Z1
    if (!isBinary()) {
      Z1 = G - X * (X.transpose() * X).ldlt().solve(X.transpose() * G);
    } else {
      Eigen::VectorXf v_sqrt = v.cwiseSqrt();
      Z1 = v_sqrt.asDiagonal() * G - v_sqrt.asDiagonal() * X * ( X.transpose() * v_sqrt.asDiagonal() * X ).ldlt().solve(X.transpose() * v_sqrt.asDiagonal() * G);
    }
    Z1 = Z1 /sqrt(2); // follow SKAT R package convention (divides sqrt{2})


    // calculate labmda.rho
    lambdas.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      Eigen::LDLT<Eigen::MatrixXf> chol;
      chol.compute(R_rhos[i]); // R_rhos[i] = L * L^T
      Eigen::MatrixXf Z2 = Z1 * chol.matrixL();
      Eigen::MatrixXf K = Z2.transpose() * Z2;
      getEigen(K, &lambdas[i]);

      if (i == 1) {
        Eigen::MatrixXf L = chol.matrixL();
        dumpEigen("R_rhos", R_rhos[i]);
        dumpEigen("chol.L", L);
        dumpEigen("Z2", Z2);        
        dumpEigen("K", K);
      }
    }

    // calculate some parameters (for Z(I-M)Z part)
    Eigen::MatrixXf z_bar = Z1.rowwise().sum() / Z1.cols();
    double z_norm = z_bar.squaredNorm();
    Eigen::MatrixXf ZMZ = z_bar.transpose() * Z1;
    ZMZ = (ZMZ.transpose() * ZMZ) / z_norm;
    Eigen::MatrixXf ZIMZ = Z1.transpose() * Z1 - ZMZ; // Z(I-M)Z' = ZZ' - ZMZ'
    getEigen(ZIMZ, &this->lambda);
    this->VarZeta  = 4.0 * (ZMZ.array() * ZIMZ.array()).sum();
    this->MuQ = lambda.sum();
    this->VarQ = 2.0 * (lambda.array() * lambda.array()).sum()  + VarZeta;

    double temp = (lambda.array() * lambda.array()) .sum() ;
    double KerQ = (lambda.array() * lambda.array() * lambda.array() * lambda.array()).sum() /  temp / temp * 12;
    this->Df = 12/ KerQ;

    // calculate tau
    taus.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      taus[i] = nMarker * nMarker * rhos[i] * z_norm + (1.0 - rhos[i]) * (z_bar.transpose() * Z1).array().square().sum() / z_norm;
    }

    // calculate moments
    moments.resize(nRho);
    for (int i = 0; i < nRho; ++i) {
      getMoment(lambdas[i], &moments[i]);
    }

    // calculate min p
    double minP;
    int minIndex;
    for (int i = 0; i < nRho; ++i) {
      if (i == 0) {
        minP = getPvalByMoment(Qs[i], moments[i]);
        minIndex = i;
      } else {
        double p = getPvalByMoment(Qs[i], moments[i]);
        if (p < minP) {
          minP = p;
          minIndex = i;
        }
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
    gsl_function F;
    F.function = integrandDavies;
    F.params = this;
    if (integration.integrate(F)) {
      fprintf(stderr, "%s:%d integration failed\n", __FILE__, __LINE__);
    };
    this->pValue = 1.0 - integration.getResult();

    return 0;
  }

  double getPvalDavies(double Q, const Eigen::MatrixXf& lambda) {
    double pval = 0.0;
    this->mixChiSq.reset();
    for (int i = 0; i < lambda.rows(); ++i) {
      this->mixChiSq.addLambda(lambda(i, 0));
    }
    pval = this->mixChiSq.getPvalue(Q);
    return pval;
  }

  double getPvalLiu(double Q, const Eigen::MatrixXf& lambda) {
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
      double v = (Qs_minP[i] - taus[i]) / (1.0 - rhos[i]);
      if (v < kappa) {
        kappa = v;
      }
    }
    double temp;
    if (kappa > lambda.sum() * 10000) {
      temp = 0.0;
    } else {
      double Q = (kappa - MuQ) * sqrt (VarQ - VarZeta) / sqrt(VarQ) + MuQ;
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
      double v = (Qs_minP[i] - taus[i]) / (1.0 - rhos[i]);
      if (v < kappa) {
        kappa = v;
      }
    }
    double Q = (kappa - MuQ) / sqrt(VarQ) * sqrt(2.0 * Df) + Df;
    return (gsl_cdf_chisq_P(Q, Df) * gsl_ran_chisq_pdf(x, 1.0));
  }
  double GetPvalue() const { return this->pValue; };
  double GetQ() const { return this->Q; };
  double GetRho() const {return this->rho;}

 private:
  bool isBinary() {
    return false;
  }
  int getRrho(double rho, Eigen::MatrixXf* R, int dim) {
    (*R).setConstant(dim, dim, rho);
    for (int i = 0; i < dim; ++i) {
      (*R)(i,i) = 1.0;
    }
    return 0;
  }
  int getEigen(Eigen::MatrixXf& k, Eigen::MatrixXf* lambda) {
    es.compute(k);
    Eigen::VectorXf values = es.eigenvalues();
    int n = values.size();
    int numNonZero = 0;
    float sumNonZero = 0.;
    for (int i = 0; i < n; ++i) {
      if (values(i) > 0) {
        ++ numNonZero ;
        sumNonZero += values(i);
      }
    }
    if (numNonZero == 0) {
      return -1;
    }
    float t = sumNonZero / numNonZero / 100000;
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
      (*lambda)(i , 0) =  values(n - 1 - i) ;
    }
    
    return 0;
  }
  int getMoment(const Eigen::MatrixXf& lambda, Moment* moment) {
    Moment& m = *moment;
    double c[4];
    Eigen::ArrayXf la = lambda.array();
    c[0] = la.sum();
    c[1] = (la * la).sum();
    c[2] = (la * la * la).sum();
    c[3] = (la * la* la* la).sum();

    m.muQ = c[0];
    double sigmaQ = sqrt(2 *c[1]);
    double s1 = c[2] / c[1] / sqrt(c[1]);
    double s2 = c[3] / (c[1] * c[1]);

    // double beta1 = sqrt(8.)*s1;
    // double beta2 = 12*s2;

    double a, d, l;
    if(s1 * s1 > s2){
      a = 1/(s1 - sqrt(s1*s1 - s2));
      d = (s1 * a - 1.0 *a*a);
      l = a*a - 2*d;
    } else {
      // type1 = 1;
      l = 1./s2;
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
    double pval = gsl_cdf_chisq_Q(Q_Norm,  df);
    return pval;
  }

  double getQvalByMoment(double min_pval, const Moment& moment) {
    double  muQ  = moment.muQ;
    double  varQ = moment.varQ;
    double  df   = moment.df;
    double  q_org = gsl_cdf_chisq_Qinv(min_pval, df);
    double  q = (q_org - df) / sqrt(2. * df) * sqrt(varQ) + muQ;
    return q;
  }
 private:
  Eigen::VectorXf res;     // residual
  Eigen::VectorXf v;
  Eigen::MatrixXf X;
  Eigen::MatrixXf G;
  Eigen::VectorXf w;     // residual

  Eigen::MatrixXf Z1;
  std::vector< float> rhos;
  std::vector< Eigen::MatrixXf > R_rhos;
  std::vector< Eigen::MatrixXf > lambdas;
  std::vector<Moment> moments;
  std::vector<float> Qs;
  std::vector<float> Qs_minP;
  std::vector<float> taus;

  int nPeople;
  int nMarker;
  int nCovariate;
  int nRho;

  MixtureChiSquare mixChiSq;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;

  double MuQ;
  double VarQ;
  double VarZeta;
  double Df;
  Eigen::MatrixXf lambda;

  double pValue;
  double Q;
  double rho;
};
SkatO::SkatO() { this->skatoImpl = new SkatOImpl; }
SkatO::~SkatO() { delete this->skatoImpl; }
void SkatO::Reset() { this->skatoImpl->Reset(); }

int SkatO::Fit(
    Vector& res_G,  // residual under NULL -- may change when permuting
    Vector& v_G,    // variance under NULL -- may change when permuting
    Matrix& X_G,    // covariance
    Matrix& G_G,    // genotype
    Vector& w_G)    // weight
{
  return this->skatoImpl->Fit(res_G, v_G, X_G, G_G, w_G);
};

double SkatO::GetPvalue() const { return this->skatoImpl->GetPvalue(); }

double SkatO::GetQ() const { return this->skatoImpl->GetQ(); }

double SkatO::GetRho() const { return this->skatoImpl->GetRho(); }

double integrandDavies(double x, void* param) {
  SkatO::SkatOImpl* p = (SkatO::SkatOImpl*) param;
  return p->computeIntegrandDavies(x);
}
double integrandLiu(double x, void* param) {
  SkatO::SkatOImpl* p = (SkatO::SkatOImpl*) param;
  return p->computeIntegrandLiu(x);
}
