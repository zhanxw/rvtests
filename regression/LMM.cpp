#include "LMM.h"

#include <limits>
#include <vector>

#include "nlopt.h"

#include "Eigen/Dense"
#include "EigenMatrix.h"
#include "EigenMatrixInterface.h"

#include "MathMatrix.h"
#define PI 3.141592653

double goalFunction(unsigned n, const double* x, double* grad,
                    void* my_func_data);

#undef DEBUG
// #define DEBUG

class LMM::Impl {
 public:
  Impl() : nSample(-1), U(-1), V(-1), llk(NAN) {}

  int FitNullModel(Matrix& Xnull, Matrix& Y) {
    G_to_Eigen(Xnull, &this->x);
    G_to_Eigen(Y, &this->y);

    // append identity matrix
    AppendIdentity();

    // initialize beta and sigma2
    sigma2.resize(kinship.size());
    for (size_t i = 0; i < sigma2.size(); ++i) {
      sigma2[i] = 1.0;
    }
    calculateSigmaMat();
    calculateBeta();

    // fit null model
    calculateLLK();
    double oldLLK = llk;
    int time = 1;
    double diff;
    while (true) {
      calculateSigma2();
      diff = llk - oldLLK;
      if (diff < 1e-6) {
#ifdef DEBUG
        fprintf(stderr, "Model converges or llk cannot be improved\n");
#endif
        break;
      }
      oldLLK = llk;
      time++;
      if (time > 100) {
#ifdef DEBUG
        fprintf(stderr, "Model probably do not converge at 100 times, llk = %g",
                llk);
#endif
        break;
      }
    }
    return 0;
  }
  int TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol) {
    G_to_Eigen(Xcol, &this->g);
    U = (resid.transpose() * sigmaMatInv * g).sum();
    V_GG = (g.transpose() * sigmaMatInv * g);
    V_XG = (x.transpose() * sigmaMatInv * g);
    V_XX = (x.transpose() * sigmaMatInv * x);
    const int m = V_XX.rows();
    V = (V_GG -
         V_XG.transpose() * V_XX.ldlt().solve(Eigen::MatrixXf::Identity(m, m)) *
             V_XG)
            .sum();

    return 0;
  }
  int AppendIdentity() {
    int n = x.rows();
    AppendKinship(Eigen::MatrixXf::Identity(n, n));
    return 0;
  }
  int AppendKinship(const Eigen::MatrixXf& k) {
    this->kinship.push_back(k);
    return 0;
  }
  int AppendKinship(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
    Eigen::MatrixXf mat =
        kinshipU.mat * kinshipS.mat * kinshipU.mat.transpose();
    return AppendKinship(mat);
  }
  int AppendKinship(const EigenMatrix& k) { return AppendKinship(k.mat); }
  int AppendKinship(Matrix& k) {
    Eigen::MatrixXf mat;
    G_to_Eigen(k, &mat);
    return AppendKinship(mat);
  }
  // for Score Test
  double GetUStat() const { return this->U; }
  double GetVStat() const { return this->V; }
  const std::vector<double>& GetSigma2() { return this->sigma2; }

  void setSigma2(unsigned n, const double* x) {
    this->sigma2.resize(n);
    for (unsigned i = 0; i < n; ++i) {
      this->sigma2[i] = x[i];
    }
  }
  void calculateSigmaMat() {
    if (kinship.empty()) return;
    nSample = kinship[0].rows();
    sigmaMat.noalias() = sigma2[0] * (kinship[0]);
    for (size_t i = 1; i < kinship.size(); ++i) {
      sigmaMat.noalias() += sigma2[i] * (kinship[i]);
    }
    cholesky = sigmaMat.llt();
    sigmaMatInv = cholesky.solve(Eigen::MatrixXf::Identity(nSample, nSample));
  }
  void calculateBeta() {
    beta = (x.transpose() * sigmaMatInv * x)
               .ldlt()
               .solve(x.transpose() * sigmaMatInv * y);
    resid = y - x * beta;
  }
  void calculateSigma2() {
    // use non linear optimization to calculate sigma2
    nlopt_opt opt;
    const int nParam = sigma2.size();
    std::vector<double> lb(nParam, 1e-10);

    opt = nlopt_create(NLOPT_LN_COBYLA, nParam);
    nlopt_set_lower_bounds(opt, lb.data());
    nlopt_set_min_objective(opt, goalFunction, this);
    nlopt_set_xtol_rel(opt, 1e-4);

    double minf; /* the minimum objective value, upon return */
    int retCode = nlopt_optimize(opt, sigma2.data(), &minf);
    if (retCode < 0) {
#ifdef DEBUG
      printf("nlopt failed [ %d ]!\n", retCode);
#endif
    } else {
#ifdef DEBUG
      printf("found minimum at %g\n", minf);
#endif
    }
    nlopt_destroy(opt);
  }
  void calculateLLK() {
    llk = (nSample * log(2.0 * PI));
    const Eigen::MatrixXf& choleskyMatrixL = cholesky.matrixL();
    llk += (choleskyMatrixL.diagonal().array().log().sum()) * 2.0;
    llk += (resid.transpose() * sigmaMatInv * resid).sum();
    llk *= -0.5;
  }
  double getLLK() const { return this->llk; }

 private:
  // Vector sigma2_g;
  std::vector<double> sigma2;
  int nSample;
  double U;
  double V;
  double llk;
  Eigen::MatrixXf sigmaMat;
  Eigen::MatrixXf sigmaMatInv;
  Eigen::MatrixXf x;  // covariate
  Eigen::VectorXf beta;
  Eigen::VectorXf resid;
  Eigen::MatrixXf y;
  Eigen::MatrixXf g;
  Eigen::MatrixXf V_XX;
  Eigen::MatrixXf V_GG;
  Eigen::MatrixXf V_XG;
  std::vector<Eigen::MatrixXf> kinship;
  Eigen::LLT<Eigen::MatrixXf> cholesky;
};

// return minus llk
double goalFunction(unsigned n, const double* x, double* grad,
                    void* my_func_data) {
  LMM::Impl* p = (LMM::Impl*)my_func_data;

#ifdef DEBUG
  fprintf(stderr, "In goalFunction()\n");
// fprintf(stderr, "set sigma2\n");
#endif
  // p->setSigma2(n, x);
  // fprintf(stderr, "calculate SimgaMat\n");
  p->calculateSigmaMat();
  // fprintf(stderr, "caluculate beta\n");
  p->calculateBeta();
  // fprintf(stderr, "calculate LLK\n");
  p->calculateLLK();
// fprintf(stderr, "return llk\n");
#ifdef DEBUG
  fprintf(stderr, "sigma2 = ");
  for (unsigned i = 0; i != n; ++i) {
    fprintf(stderr, "%g ", x[i]);
  }
  fprintf(stderr, " llk = %g\n", -p->getLLK());
#endif
  return -p->getLLK();
}

//////////////////////////////////////////////////
// LMM Interface
//////////////////////////////////////////////////
LMM::LMM() { this->impl = new Impl; }
LMM::~LMM() { delete this->impl; }

int LMM::FitNullModel(Matrix& Xnull, Matrix& y) {
  return this->impl->FitNullModel(Xnull, y);
}
int LMM::TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol) {
  return this->impl->TestCovariate(Xnull, y, Xcol);
}
int LMM::AppendKinship(const EigenMatrix& kinshipU,
                       const EigenMatrix& kinshipS) {
  return this->impl->AppendKinship(kinshipU, kinshipS);
}
int LMM::AppendKinship(const EigenMatrix& kinship) {
  return this->impl->AppendKinship(kinship);
}
int LMM::AppendKinship(Matrix& kinship) {
  return this->impl->AppendKinship(kinship);
}

// for Score Test
double LMM::GetUStat() const { return this->impl->GetUStat(); }
double LMM::GetVStat() const { return this->impl->GetVStat(); }
const std::vector<double>& LMM::GetSigma2() const {
  return this->impl->GetSigma2();
}
