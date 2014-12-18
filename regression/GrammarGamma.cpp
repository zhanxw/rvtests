#include "GrammarGamma.h"
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

#define PI 3.1415926535897

static double goalFunction(double x, void* param);

class GrammarGamma::Impl{
 public:
  Impl(){
  }
  int FitNullModel(Matrix& mat_Xnull, Matrix& mat_y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
    // type conversion
    Eigen::MatrixXf x;
    Eigen::MatrixXf y;
    G_to_Eigen(mat_Xnull, &x);
    G_to_Eigen(mat_y, &y);
    this->lambda = kinshipS.mat;
    const Eigen::MatrixXf& U = kinshipU.mat;
    // rotate
    this->ux = U.transpose() * x;
    this->uy = U.transpose() * y;

    // get beta, sigma2_g and delta
    // where delta = sigma2_e / sigma2_g
    double loglik[101];
    int maxIndex = -1;
    double maxLogLik = 0;
    for (int i = 0; i <= 100; ++i ){
      double d = exp(-10 + i * 0.2);
      getBetaSigma2(d);
      loglik[i] = getLogLikelihood(d);
      // fprintf(stderr, "%d\tdelta=%g\tll=%lf\n", i, delta, loglik[i]);
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
        // fprintf(stderr, "Minimization failed, fall back to initial guess.\n");
        this->delta = start;
      } else {
        this->delta = minimizer.getX();
        // fprintf(stderr, "minimization succeed when delta = %g, sigma2_g = %g\n", this->delta, this->sigma2_g);
      }
    }
    // store some intermediate results
    // fprintf(stderr, "maxIndex = %d, delta = %g, Try brent\n", maxIndex, delta);
    // fprintf(stderr, "beta[%d][%d] = %g\n", (int)beta.rows(), (int)beta.cols(), beta(0,0));
    this->h2 =  1.0 /(1.0 + this->delta);
    this->sigma2 = this->sigma2_g * this->h2;
    
    // we derive different formular to replace original eqn (7)
    this->gamma = (this->lambda.array() / (this->lambda.array() + this->delta)).sum() / this->sigma2_g / (this->ux.rows() - 1 ) ;
    // fprintf(stderr, "gamma = %g\n", this->gamma);
    // transformedY = \Sigma^{-1} * (y_tilda) and y_tilda = y - X * \beta
    // since \Sigma = (\sigma^2_g * h^2 ) * (U * (\lambda + delta) * U')
    // transformedY = 1 / (\sigma^2_g * h^2 ) * (U * (\lambda+delta)^{-1} * U' * (y_tilda))
    //              = 1 / (\sigma^2_g * h^2 ) * (U * \lambda^{-1} * (uResid))
    // since h^2 = 1 / (1+delta)
    // transformedY = (1 + delta/ (\sigma^2_g ) * (U * \lambda^{-1} * (uResid))
    Eigen::MatrixXf resid = y - x * (x.transpose() * x).eval().ldlt().solve(x.transpose() * y); // this is y_tilda
            
    this->transformedY.noalias() =  U.transpose() * resid;
    this->transformedY = (this->lambda.array() + this->delta).inverse().matrix().asDiagonal() * this->transformedY;
    this->transformedY = U * this->transformedY;
    this->transformedY /= this->sigma2_g;
    // fprintf(stderr, "transformedY(0,0) = %g\n", transformedY(0,0));
    
    this->ySigmaY= (resid.array() * transformedY.array()).sum();
    return 0;
  }
  int TestCovariate(Matrix& Xnull, Matrix& Y, Matrix& Xcol,
                    const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
    Eigen::MatrixXf g;
    G_to_Eigen(Xcol, &g);

    // store U'*G for computing AF later.
    const Eigen::MatrixXf& U = kinshipU.mat;
    this->ug = U.transpose() * g;

    Eigen::RowVectorXf g_mean = g.colwise().mean();
    g = g.rowwise() - g_mean;

    double gTg = g.array().square().sum();
    double t_new = (g.array() * this->transformedY.array()).sum();
    t_new = t_new * t_new / gTg;
    double t_score = t_new / this->gamma;
    this->betaG = (g.transpose() * this->transformedY).sum() / gTg / this->gamma;
    this->betaGVar = this->ySigmaY / gTg / this->gamma;

    this->pvalue = gsl_cdf_chisq_Q(t_score, 1.0);
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
    // if ( model == GrammarGamma::MLE) {
    this->sigma2_g = sumResidual2 / x.rows();
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
    // if (this->model == GrammarGamma::MLE) {
    ret = 1.0 * this->ux.rows() * log( 2.0 * PI);
    ret += (this->lambda.array() + delta).abs().log().sum();
    ret += this->ux.rows();
    ret += 1.0 * this->ux.rows() * log(this->getSumResidual2(delta));
    // }
    ret *= -0.5;
    // if (this->model == GrammarGamma::REML) {
    //   // three additional terms, however, we omit log|X'X|
    //   ret += 0.5 * this->ux.cols() * log (2.0 * PI * this->sigma2);
    //   ret -= 0.5 * log (( this->ux.transpose() * (this->lambda.array() + delta).inverse().matrix().asDiagonal() * this->ux ).determinant());
    // }
    // printf("ll = %g\n", ret);
    // this->nullLikelihood = ret;
    return ret;
  }
  // NOTE: need to fit null model fit before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) const{
    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf u1 = Eigen::MatrixXf::Ones(U.rows(), 1);
    u1 = U.transpose() *  u1;
    Eigen::MatrixXf x = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * u1;
    Eigen::MatrixXf y = (this->lambda.array() + delta).sqrt().matrix().asDiagonal() * this->ug;
    Eigen::MatrixXf beta = (x.transpose() * x).inverse() * x.transpose() * y;
    // here x is represented as 0, 1, 2, so beta(0, 0) is the mean genotype
    // multiply by 0.5 to get AF
    double af = 0.5 * beta(0, 0);
    return af;
  }
  double GetPvalue() const{
    return this->pvalue;
  }
  double GetBeta() const {
    return this->betaG;
  }
  double GetBetaVar() const {
    return this->betaGVar;
  }
  double GetSigmaG2() {
    return this->sigma2;
  }
  double GetSigmaE2() {
    return this->sigma2 * this->delta;
  }
  double GetDelta()  {
    return this->delta;
  }    // delta = sigma2_e / sigma2_g
  
   private:
  // Eigen::MatrixXf S;
  float sigma2_g;     // sigma2_g
  float delta;      // delta =  sigma2_e / sigma2_g
  float sigma2;
  float h2;
  Eigen::MatrixXf beta;
  // temporary values
  Eigen::MatrixXf uy; // U' * y
  Eigen::MatrixXf ux; // U' * x
  Eigen::MatrixXf ug; // U' * g
  Eigen::MatrixXf lambda;
  // for score test
  // Eigen::MatrixXf uResid; // U' * (y - mean(y))
  double gamma;
  Eigen::MatrixXf transformedY;
  double ySigmaY;
  double pvalue;
  double betaG;
  double betaGVar;
  double sumResidual2; // sum (  (Uy - Ux *beta)^2/(lambda + delta) )
};

//////////////////////////////////////////////////
// GrammarGamma Interface
//////////////////////////////////////////////////
GrammarGamma::GrammarGamma(){
  this->impl = new Impl();
}
GrammarGamma::~GrammarGamma(){
  delete this->impl;
}

// @return 0 when success
int GrammarGamma::FitNullModel(Matrix& Xnull, Matrix& y,
                               const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
  return this->impl->FitNullModel(Xnull, y, kinshipU, kinshipS);
}

int GrammarGamma::TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol,
                                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
  return this->impl->TestCovariate(Xnull, y, Xcol, kinshipU, kinshipS);
}
double GrammarGamma::GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS){
  return this->impl->GetAF(kinshipU, kinshipS);
}
double GrammarGamma::GetPvalue(){
  return this->impl->GetPvalue();
}
double GrammarGamma::GetBeta() { return this->impl->GetBeta();};
double GrammarGamma::GetBetaVar() { return this->impl->GetBetaVar();};
double GrammarGamma::GetSigmaE2() const { return this->impl->GetSigmaE2(); };
double GrammarGamma::GetSigmaG2() const { return this->impl->GetSigmaG2(); };
double GrammarGamma::GetDelta()  const { return this->impl->GetDelta(); };    // delta = sigma2_e / sigma2_g

// need to negaive the MLE to minize it
double goalFunction(double x, void* param) {
  GrammarGamma::Impl* p = (GrammarGamma::Impl*) param;
  p->getBetaSigma2(x);
  return (- p->getLogLikelihood(x) );
}
