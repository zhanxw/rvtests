#include "FastLMM.h"
#include "EigenMatrix.h"
#include "EigenMatrixInterface.h"
#include "Eigen/Dense"
#include "GSLMinimizer.h"
#include "gsl/gsl_cdf.h" // use gsl_cdf_chisq_Q

#if 0
#include <fstream>
void dumpToFile(const Eigen::MatrixXf& mat, const char* fn) {
  std::ofstream out(fn);
  out << mat.rows() << "\t" << mat.cols() << "\n";
  out << mat;
  out.close();
}
#endif
// #define EIGEN_NO_DEBUG
#undef DEBUG
// #define DEBUG

#define PI 3.1415926535897

static double goalFunction(double x, void* param);

class FastLMM::Impl{
 public:
  Impl(FastLMM::Test test, FastLMM::Model model): test(test), model(model){
  }
  int FitNullModel(Matrix& mat_Xnull, Matrix& mat_y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
    // sanity check
    if (mat_Xnull.rows != mat_y.rows) return -1;
    if (mat_Xnull.rows != kinshipU.mat.rows()) return -1;
    if (mat_Xnull.rows != kinshipS.mat.rows()) return -1;
    
    // type conversion
    G_to_Eigen(mat_Xnull, &this->ux);
    G_to_Eigen(mat_y, &this->uy);
    this->lambda = kinshipS.mat;
    const Eigen::MatrixXf& U = kinshipU.mat;
    // rotate
    this->ux = U.transpose() * this->ux;
    this->uy = U.transpose() * this->uy;

    // get beta, sigma and delta
    // where delta = sigma2_e / sigma2_g
    double loglik[101];
    int maxIndex = -1;
    double maxLogLik = 0;
    for (int i = 0; i <= 100; ++i ){
      delta = exp(-10 + i * 0.2);
      getBetaSigma2(delta);
      loglik[i] = getLogLikelihood(delta);
#ifdef DEBUG
      fprintf(stderr, "%d\tdelta=%g\tll=%lf\n", i, delta, loglik[i]);
      fprintf(stderr, "beta(0)=%lf\tsigma2=%lf\n",
              beta(0), sigma2);
#endif
      if (std::isnan(loglik[i])) {
        continue;
      }
      if (maxIndex < 0 || loglik[i] > maxLogLik) {
        maxIndex = i;
        maxLogLik = loglik[i];
      }
    }
    if (maxIndex < -1) {
      fprintf(stderr, "Cannot optimize\n");
      return -1;
    }
#ifdef DEBUG
    fprintf(stderr, "maxIndex = %d\tll=%lf\t\tbeta(0)=%lf\tsigma2=%lf\n",
            maxIndex, maxLogLik, beta(0), sigma2);
#endif

    if (maxIndex == 0 || maxIndex == 100) {
      // on the boundary
      // do not try maximize it.
    } else {
      gsl_function F;
      F.function = goalFunction;
      F.params = this;

      Minimizer minimizer;
      double lb = exp(-10 + (maxIndex-1) * 0.2);
      double ub = exp(-10 + (maxIndex+1) * 0.2);
      double start =  exp(-10 + maxIndex * 0.2);
      if (minimizer.minimize(F, start, lb, ub)) {
        fprintf(stderr, "Minimization failed, fall back to initial guess.\n");
        this->delta = start;
      } else {
        this->delta = minimizer.getX();
#ifdef DEBUG       
        fprintf(stderr, "minimization succeed when delta = %g, sigma2 = %g\n", this->delta, this->sigma2);
#endif
      }
    }
    // store some intermediate results
#ifdef DEBUG       
    fprintf(stderr, "maxIndex = %d, delta = %g, Try brent\n", maxIndex, delta);
    fprintf(stderr, "beta[%d][%d] = %g\n", (int)beta.rows(), (int)beta.cols(), beta(0,0));
#endif
    if (this->test == FastLMM::LRT) {
      this->nullLikelihood = getLogLikelihood(this->delta);
    } else if (this->test == FastLMM::SCORE) {
      this->uResid = this->uy - this->ux * this->beta;
    }
    return 0;
  }
  int TestCovariate(Matrix& Xnull, Matrix& Y, Matrix& Xcol,
                    const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
    // obtain U
    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf g;
    G_to_Eigen(Xcol, &g);
    this->ug = U.transpose() * g;

    // depends on LRT test or Score tests
    if (this->test == FastLMM::LRT) {
      // need to fit alternative model
      // 1. assign obs and parameters

      Eigen::MatrixXf altUx; // X under alternative model
      // cbind_G_to_Eigen(Xnull, Xcol, &altUx);
      altUx.resize(this->ux.rows(), this->ux.cols() + this->ug.cols());
      altUx << this->ux, this->ug;
      
      Eigen::MatrixXf& altUy =  this->uy;
      Eigen::MatrixXf& altLambda =  this->lambda;
      // 2. estimate beta and sigma2 using delta under the null model
      Eigen::MatrixXf x = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * altUx;
      Eigen::MatrixXf y = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * altUy;
      Eigen::MatrixXf altBeta = (x.transpose() * x).eval().ldlt().solve(x.transpose() * y);
      double altSumResidual2 = (( altUy.array() - (altUx * altBeta).array() ).square() / (altLambda.array() + delta)).sum();

      double altSigma2;
      if ( model == FastLMM::MLE) {
        altSigma2 = altSumResidual2 / x.rows();
      } else {
        altSigma2 = altSumResidual2 / (x.rows() - x.cols());
      }
      // 3. get new likelhood
      double ret = 0;
      if (this->model == FastLMM::MLE) {
        ret = 1.0 * altUx.rows() * log( 2.0 * PI);
        ret += (altLambda.array() + delta).log().sum();
        ret += altUx.rows();
        ret += 1.0 * altUx.rows() * log(altSumResidual2);
      }
      ret *= -0.5;
      if (this->model == FastLMM::REML) {
        // three additional terms, however, we omit log|X'X|
        ret += 0.5 * altUx.cols() * log (2.0 * PI * altSigma2);
        ret -= 0.5 * log (( altUx.transpose() * (altLambda.array() + delta).inverse().matrix().asDiagonal() * altUx ).determinant());
      }
      this->altLikelihood = ret;
      this->stat = 2.0 * (this->altLikelihood - this->nullLikelihood); // stat ~ X^2 1df distribution
      this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);
      return 0;
    }

    if (this->test == FastLMM::SCORE) {
      // just return score test statistics
      // Eigen::RowVectorXf g_mean = g.colwise().mean();
      Eigen::MatrixXf u_g_center = U.transpose() *  (g.rowwise() - g_mean);
      this->Ustat = (  (  (u_g_center).array() *
                          (  (this->uResid)).array() ) /
                       (lambda.array() + delta)
                       ).sum() / this->sigma2;
      this->Vstat = ((u_g_center).array().square() / (lambda.array() + delta)).sum() / this->sigma2 ;
      if (this->Vstat > 0.0) {
        this->stat = this->Ustat * this->Ustat / this->Vstat;
        this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);
      } else {
        this->stat = 0;
        this->pvalue = 1.0;
      }
      return 0;
    }
    return 0;
  }
  double getSumResidual2(double delta) {
    return (( this->uy.array() - (this->ux * this->beta).array() ).square() / (this->lambda.array() + delta)).sum();
  }
  void getBetaSigma2(double delta) {
    Eigen::MatrixXf x = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * this->ux;
    Eigen::MatrixXf y = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * this->uy;
    this->beta = (x.transpose() * x).eval().ldlt().solve(x.transpose() * y);
    double sumResidual2 = getSumResidual2(delta);
    if ( model == FastLMM::MLE) {
      this->sigma2 = sumResidual2 / x.rows();
    } else {
      this->sigma2 = sumResidual2 / (x.rows() - x.cols());
    }
  }
  /**
   * NOTE: it's necesary to calculate beta and sigma2 beforehand.
   * use estimated beta(delta) and estimated sigma2(delta)
   * to calculate the likelihood of given delta
   */
  double getLogLikelihood(double delta) {
    double ret = 0;
    if (this->model == FastLMM::MLE) {
      ret = 1.0 * this->ux.rows() * log( 2.0 * PI);
      ret += (this->lambda.array() + delta).abs().log().sum();
      ret += this->ux.rows();
      ret += 1.0 * this->ux.rows() * log(this->getSumResidual2(delta));
    }
    ret *= -0.5;
    if (this->model == FastLMM::REML) {
      // three additional terms, however, we omit log|X'X|
      ret += 0.5 * this->ux.cols() * log (2.0 * PI * this->sigma2);
      ret -= 0.5 * log (( this->ux.transpose() * (this->lambda.array() + delta).inverse().matrix().asDiagonal() * this->ux ).determinant());
    }
    // printf("ll = %g\n", ret);
    this->nullLikelihood = ret;
    return ret;
  }
  // NOTE: need to fit null model before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) const{
    // K = U S U^t = U * (lambda + delta) * U^t
    // 1 = (n by 1)  matrix
    // G = (n by 1) genotype matrix
    // AF = (1^t K^(-1) 1)^(-1) * (1^t K^(-1) G)
    //    = (1^t U (lambda + delta)^(-1) * U^t * 1)^(-1) *
    //      (1^t U (lambda + delta)^(-1) * U^t * G)
    //    = ((U^t*1)^t * (lambda + delta)^(-1) * (U^t*1))^(-1)
    //      ((U^t*1)^t * (lambda + delta)^(-1) * (U^t * G))
    //    = (u1s * u1)^(-1) * (u1s * ug)
    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf u1 = U.transpose().rowwise().sum();
    // Eigen::MatrixXf ug = U.transpose() .rowwise().sum();
    // Eigen::MatrixXf u1s = u1.transpose() * (this->lambda.array() + delta).matrix().asDiagonal();
    Eigen::ArrayXf u1s = u1.array() / (this->lambda.array() + delta).array();
    double denom = (u1s * u1.array()).sum();
    if (denom == 0.0) {
      return 0.0;
    }
    double numer = (u1s * ug.array()).sum();
    double beta =  numer / denom;

    // here x is represented as 0, 1, 2, so beta(0, 0) is the mean genotype
    // multiply by 0.5 to get AF
    double af = 0.5 * beta;
    return af;
  }
  // NOTE: need to fit null model before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS, Matrix& Xcol) const{
    Eigen::MatrixXf g;
    G_to_Eigen(Xcol, &g);

    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf u1 = U.transpose().rowwise().sum();
    Eigen::MatrixXf ug = U.transpose() * g;
    Eigen::ArrayXf u1s = u1.array() / (this->lambda.array() + delta).abs().array();
    // This is usually the same, must avoid future calculations...
    double denom = (u1s * u1.array()).sum();
    if (denom == 0.0) {
      return 0.0;
    }
    double numer = (u1s * ug.array()).sum();
    double beta =  numer / denom;

    // here x is represented as 0, 1, 2, so beta(0, 0) is the mean genotype
    // multiply by 0.5 to get AF
    double af = 0.5 * beta;
    return af;
  }
  // NOTE: need to fit null model before calling this function
  // NOTE2: assuming kinship matrices are unchanged
  double FastGetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS, Matrix& Xcol) const{
    const Eigen::MatrixXf& U = kinshipU.mat;
    static initialized = false;
    static Eigen::MatrixXf g;
    static Eigen::MatrixXf u1s;
    static double denom;

    if (!initialized) {
      Eigen::MatrixXf u1 = U.transpose().rowwise().sum();
      Eigen::MatrixXf ug = U.transpose() * g;
      u1s = u1.array() / (this->lambda.array() + delta).abs().array();
      denom = (u1s * u1.array()).sum();
      initialized = true;
    }
    
    if (denom == 0.0) {
      return 0.0;
    }
    G_to_Eigen(Xcol, &g);
    double numer = (u1s * ug.array()).sum();
    double beta =  numer / denom;

    // here x is represented as 0, 1, 2, so beta(0, 0) is the mean genotype
    // multiply by 0.5 to get AF
    double af = 0.5 * beta;
    return af;
  }
  
  double GetPvalue() const{
    return this->pvalue;
  }
  double GetUStat() const {
    return this->Ustat;
  }
  double GetVStat() const {
    return this->Vstat;
  }
  double GetEffect() {
    if (this->Vstat != 0.0)
      return this->Ustat / this->Vstat;
    return 0.;
  };   // U/V
  double GetSE()  {
    if (this->Vstat != 0.0)
      return 1 / sqrt(this->Vstat);
    return 0.;
  };   // 1/sqrt(V)
  double GetSigmaG2() {
    return this->sigma2;
  };
  double GetSigmaE2() {
    return this->sigma2 * this->delta;
  };
  double GetDelta()  {
    return this->delta;
  };    // delta = sigma2_e / sigma2_g
  
  double GetNullLogLikelihood() const {
    return this->nullLikelihood;
  }
  double GetAltLogLikelihood() const {
    return this->altLikelihood;
  }

 private:
  // Eigen::MatrixXf S;
  float sigma2;     // sigma2_g
  float delta;      // delta =  sigma2_e / sigma2_g
  Eigen::MatrixXf beta;
  // temporary values
  Eigen::MatrixXf uy; // U' * y
  Eigen::MatrixXf ux; // U' * x
  Eigen::MatrixXf ug; // U' * g
  Eigen::MatrixXf lambda;
  // for LRT test
  double nullLikelihood;
  double altLikelihood;
  // for score test
  Eigen::MatrixXf uResid; // U' * (y - mean(y))
  double Ustat;
  double Vstat;
  double stat;
  // for both tests
  double pvalue;
  //
  FastLMM::Test test;
  FastLMM::Model model;
  double sumResidual2; // sum (  (Uy - Ux *beta)^2/(lambda + delta) )
};

//////////////////////////////////////////////////
// FastLMM Interface
//////////////////////////////////////////////////
FastLMM::FastLMM(Test test, Model model){
  this->impl = new Impl(test, model);
}
FastLMM::~FastLMM(){
  delete this->impl;
}

// @return 0 when success
int FastLMM::FitNullModel(Matrix& Xnull, Matrix& y,
                          const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
  return this->impl->FitNullModel(Xnull, y, kinshipU, kinshipS);
}
int FastLMM::TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol,
                           const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
  return this->impl->TestCovariate(Xnull, y, Xcol, kinshipU, kinshipS);
}
double FastLMM::GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
  return this->impl->GetAF(kinshipU, kinshipS);
}
double FastLMM::GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS, Matrix& Xcol){
  return this->impl->GetAF(kinshipU, kinshipS, Xcol);
}
double FastLMM::GetPvalue(){
  return this->impl->GetPvalue();
}
double FastLMM::GetUStat() { return this->impl->GetUStat();};
double FastLMM::GetVStat() { return this->impl->GetVStat();};
double FastLMM::GetEffect() { return this->impl->GetEffect(); };   // U/V
double FastLMM::GetSE()      { return this->impl->GetSE(); };   // 1/sqrt(V)
double FastLMM::GetSigmaE2() { return this->impl->GetSigmaE2(); };
double FastLMM::GetSigmaG2() { return this->impl->GetSigmaG2(); };
double FastLMM::GetDelta()  { return this->impl->GetDelta(); };    // delta = sigma2_e / sigma2_g
double FastLMM::GetNullLogLikelihood() { return this->impl->GetNullLogLikelihood(); };
double FastLMM::GetAltLogLikelihood() { return this->impl->GetAltLogLikelihood(); };


// need to negaive the MLE to minize it
double goalFunction(double x, void* param) {
  FastLMM::Impl* p = (FastLMM::Impl*) param;
  p->getBetaSigma2(x);
  return (- p->getLogLikelihood(x) );
}
