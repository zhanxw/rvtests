#include "MetaCov.h"
#include "EigenMatrix.h"
#include "EigenMatrixInterface.h"
#include "Eigen/Dense"
#include "GSLMinimizer.h"

#undef DEBUG
// #define DEBUG
#define PI 3.1415926535897

static double goalFunction(double x, void* param);

class MetaCov::Impl{
 public:
  Impl(){
  }
  int FitNullModel(Matrix& mat_Xnull, Matrix& mat_y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
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
    // if (this->test == MetaCov::LRT) {
    // this->nullLikelihood = getLogLikelihood(this->delta);
    // } else if (this->test == MetaCov::SCORE) {
    //   this->uResid = this->uy - this->ux * this->beta;
    // }
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
    // if ( model == MetaCov::MLE) {
      this->sigma2 = sumResidual2 / x.rows();
    // } else {
    //   this->sigma2 = sumResidual2 / (x.rows() - x.cols());
    // }
  }
  /**
   * NOTE: it's necesary to calculate beta and sigma2 beforehand.
   * use estimated beta(delta) and estimated sigma2(delta)
   * to calculate the likelihood of given delta
   */
  double getLogLikelihood(double delta) {
    double ret = 0;
    // if (this->model == MetaCov::MLE) {
      ret = 1.0 * this->ux.rows() * log( 2.0 * PI);
      ret += (this->lambda.array() + delta).abs().log().sum();
      ret += this->ux.rows();
      ret += 1.0 * this->ux.rows() * log(this->getSumResidual2(delta));
    // }
    ret *= -0.5;
    // if (this->model == MetaCov::REML) {
    //   // three additional terms, however, we omit log|X'X|
    //   ret += 0.5 * this->ux.cols() * log (2.0 * PI * this->sigma2);
    //   ret -= 0.5 * log (( this->ux.transpose() * (this->lambda.array() + delta).inverse().matrix().asDiagonal() * this->ux ).determinant());
    // }
    // printf("ll = %g\n", ret);
    // this->nullLikelihood = ret;
    return ret;
  }
  // // NOTE: need to fit null model fit before calling this function
  // double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) const{
  //   const Eigen::MatrixXf& U = kinshipU.mat;
  //   Eigen::MatrixXf u1 = Eigen::MatrixXf::Ones(U.rows(), 1);
  //   u1 = U.transpose() *  u1;
  //   Eigen::MatrixXf x = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * u1;
  //   Eigen::MatrixXf y = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * this->ug;
  //   Eigen::MatrixXf beta = (x.transpose() * x).inverse() * x.transpose() * y;
  //   // here x is represented as 0, 1, 2, so beta(0, 0) is the mean genotype
  //   // multiply by 0.5 to get AF
  //   double af = 0.5 * beta(0, 0);
  //   return af;
  // }
  // U * ( x - center(x) )
  int TransformCentered(std::vector<double>* geno,
                        const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
    // type conversion
    int n = geno->size();
    this->g.resize(n, 1);
    for(int i = 0; i < n; ++i){
      this->g(i) = (*geno)[i];
    }
    Eigen::RowVectorXf g_mean = g.colwise().mean();
    const Eigen::MatrixXf& U = kinshipU.mat;
    this->g = U * (g.rowwise() - g_mean);
    for(int i = 0; i < n; ++i){
      (*geno)[i] = this->g(i);
    }
    return 0;
  }
  int GetWeight(Vector* out) const {
    const int n = lambda.rows();
    out->Dimension(n);
    for (int i = 0; i < n ; ++i){
      (*out)[i] = sigma2 * (lambda(i) + delta) ;
    }
    // fprintf(stderr, "sigma2 = %g, lambda(0) = %g, lambda(99) = %g, delta = %g\n", sigma2, lambda(0), lambda(99), delta);
    return 0;
  }

 private:
  // Eigen::MatrixXf S;
  float sigma2;     // sigma2_g
  float delta;      // delta =  sigma2_e / sigma2_g
  Eigen::MatrixXf beta;
  // temporary values
  Eigen::MatrixXf g;
  Eigen::MatrixXf uy; // U' * y
  Eigen::MatrixXf ux; // U' * x
  // Eigen::MatrixXf ug; // U' * g
  Eigen::MatrixXf lambda;
  // // for LRT test
  // double nullLikelihood;
  // double altLikelihood;
  // // for score test
  // Eigen::MatrixXf uResid; // U' * (y - mean(y))
  // double Ustat;
  // double Vstat;
  // double stat;
  // // for both tests
  // double pvalue;
  //
  double sumResidual2; // sum (  (Uy - Ux *beta)^2/(lambda + delta) )

};

//////////////////////////////////////////////////
// MetaCov Interface
//////////////////////////////////////////////////
MetaCov::MetaCov(){
  this->impl = new Impl();
}

MetaCov::~MetaCov(){
  delete this->impl;
}

// @return 0 when success
int MetaCov::FitNullModel(Matrix& Xnull, Matrix& y,
                          const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
  return this->impl->FitNullModel(Xnull, y, kinshipU, kinshipS);
}
  
int  MetaCov::TransformCentered(std::vector<double>* x,
                                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
  return this->impl->TransformCentered(x, kinshipU, kinshipS);
}

int  MetaCov::GetWeight(Vector* out){
  return this->impl->GetWeight(out);
}

// need to negaive the MLE to minize it
double goalFunction(double x, void* param) {
  MetaCov::Impl* p = (MetaCov::Impl*) param;
  p->getBetaSigma2(x);
  return (- p->getLogLikelihood(x) );
}
