#pragma GCC diagnostic ignored "-Wint-in-bool-context"

#include "FastLMM.h"

#include "libsrc/MathMatrix.h"
#include "third/eigen/Eigen/Dense"
#include "third/gsl/include/gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q

#include "EigenMatrix.h"
#include "EigenMatrixInterface.h"
#include "GSLMinimizer.h"

// set environment variable FASTLMM_DEBUG to enable debug information
// #define EIGEN_NO_DEBUG
// #undef DEBUG
// #define DEBUG

#define PI 3.1415926535897

static double goalFunction(double x, void* param);

class FastLMM::Impl {
 public:
  Impl(FastLMM::Test test, FastLMM::Model model)
      : test(test), model(model), needToCenterGentype(true) {
    FastLMM::Impl::showDebug = false;
  }
  int FitNullModel(Matrix& mat_Xnull, Matrix& mat_y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
    //
    if (std::getenv("FASTLMM_DEBUG")) {
      FastLMM::Impl::showDebug = true;
    } else {
      FastLMM::Impl::showDebug = false;
    }

    // sanity check
    if (mat_Xnull.rows != mat_y.rows) return -1;
    if (mat_Xnull.rows != kinshipU.mat.rows()) return -1;
    if (mat_Xnull.rows != kinshipS.mat.rows()) return -1;

    // type conversion
    G_to_Eigen(mat_Xnull, &this->ux);
    G_to_Eigen(mat_y, &this->uy);
    this->lambda = kinshipS.mat;
    // take absoluate value
    // otherwise, small negative lambda may be instable
    // e.g. lambd[i] = -0.00013, and delta = 8e-5, then (1/(lambda[i] + delta))
    //      becomes extremely large negative values, and it affects results
    this->lambda = this->lambda.cwiseAbs();

    const Eigen::MatrixXf& U = kinshipU.mat;

    // rotate
    this->ux = U.transpose() * this->ux;
    this->uy = U.transpose() * this->uy;

    // get beta, sigma and delta
    // where delta = sigma2_e / sigma2_g
    double loglik[101];
    int maxIndex = -1;
    double maxLogLik = 0;
    for (int i = 0; i <= 100; ++i) {
      delta = exp(-10. + i * 0.2);
      getBetaSigma2(delta);
      loglik[i] = getLogLikelihood(delta);
      if (FastLMM::Impl::showDebug) {
        fprintf(stderr, "%d\tdelta=%g\tll=%lf\t", i, delta, loglik[i]);
        fprintf(stderr, "beta(0)=%lf\tsigma2=%lf\n", beta(0), sigma2);
      }
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
    if (FastLMM::Impl::showDebug) {
      fprintf(stderr, "maxIndex = %d\tll=%lf\t\tbeta(0)=%lf\tsigma2=%lf\n",
              maxIndex, maxLogLik, beta(0), sigma2);
    }

    if (maxIndex == 0 || maxIndex == 100) {
      // on the boundary
      // do not try maximize it.
    } else {
      gsl_function F;
      F.function = goalFunction;
      F.params = this;

      Minimizer minimizer;
      double lb = exp(-10. + (maxIndex - 1) * 0.2);
      double ub = exp(-10. + (maxIndex + 1) * 0.2);
      double start = exp(-10. + maxIndex * 0.2);
      if (minimizer.minimize(F, start, lb, ub)) {
        fprintf(stderr, "Minimization failed, fall back to initial guess.\n");
        this->delta = start;
      } else {
        this->delta = minimizer.getX();
        if (FastLMM::Impl::showDebug) {
          fprintf(stderr, "minimization succeed when delta = %g, sigma2 = %g\n",
                  this->delta, this->sigma2);
          dumpToFile(kinshipS.mat, "S.mat");
          dumpToFile(kinshipU.mat, "U.mat");
        }
      }
    }
    // store some intermediate results
    if (FastLMM::Impl::showDebug) {
      fprintf(stderr, "delta = sigma2_e/sigma2_g, and sigma2 is sigma2_g\n");
      fprintf(stderr, "maxIndex = %d, delta = %g, Try brent\n", maxIndex,
              delta);
      fprintf(stderr, "beta[0][0] = %g\t sigma2_g = %g\tsigma2_e = %g\n",
              beta(0, 0), this->sigma2, delta * sigma2);
    }
    if (this->test == FastLMM::LRT) {
      this->nullLikelihood = getLogLikelihood(this->delta);
    } else if (this->test == FastLMM::SCORE) {
      this->uResid = this->uy - this->ux * this->beta;
    }

    // set scaledK
    // when K = U * S * U'
    // From  K^(-1) - K^(-1) * X * (X' * K^(-1) * X)^(-1) * X' * K^(-1)
    //     = U * (S^(-1) - S^(-1) * UX * ((UX)' * S^(-1) * (UX))^(-1) * (UX)' *
    //     S^(-1)) * U'
    // define scaledK = S^(-1) - S^(-1) * UX * ((UX)' * S^(-1) * (UX))^(-1) *
    // (UX)' * S^(-1)
    const Eigen::MatrixXf Sinv =
        (this->lambda.array() + delta).inverse().matrix().asDiagonal();
    this->scaledK = Sinv -
                    Sinv * ux * (ux.transpose() * Sinv * ux).inverse() *
                        ux.transpose() * Sinv;

    CalculateFactors(U);
    return 0;
  }
  int TestCovariate(Matrix& Xnull, Matrix& Y, Matrix& Xcol,
                    const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) {
    // obtain U
    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf g;
    G_to_Eigen(Xcol, &g);

    // depends on LRT test or Score tests
    if (this->test == FastLMM::LRT) {
      // need to fit alternative model
      // 1. assign obs and parameters
      this->ug.noalias() = U.transpose() * g;
      Eigen::MatrixXf altUx;  // X under alternative model
      // cbind_G_to_Eigen(Xnull, Xcol, &altUx);
      altUx.resize(this->ux.rows(), this->ux.cols() + this->ug.cols());
      altUx << this->ux, this->ug;

      Eigen::MatrixXf& altUy = this->uy;
      Eigen::MatrixXf& altLambda = this->lambda;
      // 2. estimate beta and sigma2 using delta under the null model
      Eigen::MatrixXf x = (this->lambda.array() + delta)
                              .sqrt()
                              .inverse()
                              .matrix()
                              .asDiagonal() *
                          altUx;
      Eigen::MatrixXf y = (this->lambda.array() + delta)
                              .sqrt()
                              .inverse()
                              .matrix()
                              .asDiagonal() *
                          altUy;
      Eigen::MatrixXf altBeta =
          (x.transpose() * x).eval().ldlt().solve(x.transpose() * y);
      double altSumResidual2 =
          ((altUy.array() - (altUx * altBeta).array()).square() /
           (altLambda.array() + delta))
              .sum();

      double altSigma2;
      if (model == FastLMM::MLE) {
        altSigma2 = altSumResidual2 / x.rows();
      } else {
        altSigma2 = altSumResidual2 / (x.rows() - x.cols());
      }
      // 3. get new likelhood
      // See 1.4
      double ret = 0;
      if (this->model == FastLMM::MLE) {
        ret = 1.0 * altUx.rows() * log(2.0 * PI);
        ret += (altLambda.array() + delta).abs().log().sum();
        ret += altUx.rows();
        ret += 1.0 * altUx.rows() * log(altSigma2);
      }
      ret *= -0.5;
      if (this->model == FastLMM::REML) {
        // three additional terms, however, we omit log|X'X|
        ret += 0.5 * altUx.cols() * log(2.0 * PI * altSigma2);
        ret -=
            0.5 *
            log((altUx.transpose() *
                 (altLambda.array() + delta).inverse().matrix().asDiagonal() *
                 altUx)
                    .determinant());
      }
      this->altLikelihood = ret;
      this->stat = 2.0 * (this->altLikelihood -
                          this->nullLikelihood);  // stat ~ X^2 1df distribution
      this->pvalue = gsl_cdf_chisq_Q(this->stat, 1.0);
      return 0;
    }

    if (this->test == FastLMM::SCORE) {
      // just return score test statistics
      Eigen::ArrayXf u_g_center;
      if (needToCenterGentype) {
        u_g_center =
            (U.transpose() * (g.rowwise() - g.colwise().mean())).eval().array();
      } else {
        u_g_center = (U.transpose() * g).eval().array();
      }
      if (FastLMM::Impl::showDebug) {
        dumpToFile(U, "U");
        dumpToFile(lambda, "lambda");
        dumpToFile(g, "g");
        dumpToFile(uResid, "uResid");
        dumpToFile(u_g_center, "u_g_center");
      }

      this->Ustat =
          (((u_g_center) * (this->uResid).array()) / (lambda.array() + delta))
              .sum() /
          this->sigma2;
      // when there is no covariate: this->Vstat = (u_g_center.square() /
      // (lambda.array() + delta)).sum() / this->sigma2 ;
      this->Vstat = ((u_g_center).matrix().transpose() * this->scaledK *
                     (u_g_center.matrix()))(0, 0) /
                    this->sigma2;
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
  /**
   * calculate U and V matrix as in score test. It's useful when @param has
   * multiple columns.
   * uMat = Xcol' * K^{-1} Y
   * vMat = {Xcol' * (K^{-1} - K^{-1} * Xnull' (Xnull' * K^{-1} * Xnull)^{-1} *
   * Xnull * K^{-1}) * Xcol} * sigma2
   */
  int CalculateUandV(Matrix& Xnull, Matrix& Y, Matrix& Xcol,
                     const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                     Matrix* uMat, Matrix* vMat) {
    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf u_g_center;
    Eigen::MatrixXf g;
    G_to_Eigen(Xcol, &g);

    if (needToCenterGentype) {
      u_g_center = (U.transpose() * (g.rowwise() - g.colwise().mean())).eval();
    } else {
      u_g_center = (U.transpose() * g).eval();
    }
    Eigen::MatrixXf u =
        (u_g_center.transpose()) *
        (lambda.array() + delta).inverse().matrix().asDiagonal() *
        (this->uResid) / this->sigma2;
    Eigen::MatrixXf v =
        (u_g_center.transpose() * this->scaledK * u_g_center) / this->sigma2;

    Eigen_to_G(u, uMat);
    Eigen_to_G(v, vMat);

#if 0
    dumpToFile(lambda, "lambda");
    dumpToFile(scaledK, "scaledK");
    dumpToFile(g, "g");
    dumpToFile(u_g_center, "u_g_center");
    dumpToFile(u, "u");
    dumpToFile(v, "v");
#endif
    return 0;
  }
  double getSumResidual2(double delta) {
    return ((this->uy.array() - (this->ux * this->beta).array()).square() /
            (this->lambda.array() + delta))
        .sum();
  }
  void getBetaSigma2(double delta) {
    // Eigen::MatrixXf x = (this->lambda.array() +
    // delta).sqrt().matrix().asDiagonal() * this->ux;
    // Eigen::MatrixXf y = (this->lambda.array() +
    // delta).sqrt().matrix().asDiagonal() * this->uy;
    // this->beta = (x.transpose() * x).eval().ldlt().solve(x.transpose() * y);
    if (FastLMM::Impl::showDebug) {
      dumpToFile(ux, "ux");
      dumpToFile(uy, "uy");
    }
    this->beta =
        (ux.transpose() *
         (this->lambda.array() + delta).abs().inverse().matrix().asDiagonal() *
         ux)
            .ldlt()
            .solve(ux.transpose() *
                   (this->lambda.array() + delta)
                       .abs()
                       .inverse()
                       .matrix()
                       .asDiagonal() *
                   uy);

    double sumResidual2 = getSumResidual2(delta);
    if (model == FastLMM::MLE) {
      this->sigma2 = sumResidual2 / ux.rows();
    } else {
      this->sigma2 = sumResidual2 / (ux.rows() - ux.cols());
    }
  }
  /**
   * NOTE: it's necesary to calculate beta and sigma2 beforehand.
   * use estimated beta(delta) and estimated sigma2(delta)
   * to calculate the likelihood of given delta
   */
  double getLogLikelihood(double delta) {
    double ret = 0;
    const double n = this->ux.rows();
    if (this->model == FastLMM::MLE) {
      ret = n * log(2.0 * PI);
      ret += (this->lambda.array() + delta).abs().log().sum();
      ret += n;
      ret += n * log(this->sigma2);
    }
    ret *= -0.5;
    if (this->model == FastLMM::REML) {
      // three additional terms, however, we omit log|X'X|
      ret += 0.5 * this->ux.cols() * log(2.0 * PI * this->sigma2);
      ret -=
          0.5 *
          log((this->ux.transpose() *
               (this->lambda.array() + delta).inverse().matrix().asDiagonal() *
               this->ux)
                  .determinant());
    }
    // printf("ll = %g\n", ret);
    this->nullLikelihood = ret;
    return ret;
  }
  // NOTE: need to fit null model before calling this function
  // use the internal variable this->ug to calculate
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS) const {
    return GetAFFromUg(kinshipU, kinshipS, this->ug);
  }
  // NOTE: need to fit null model before calling this function
  double GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
               Matrix& Xcol) const {
    Eigen::MatrixXf g;
    G_to_Eigen(Xcol, &g);

    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf ug = U.transpose() * g;
    return GetAFFromUg(kinshipU, kinshipS, ug);
  }
  // NOTE: need to fit null model before calling this function
  double GetAFFromUg(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                     const Eigen::MatrixXf& ug) const {
    // K = U S U^t = U * (lambda) * U^t
    // 1 = (n by 1)  matrix
    // G = (n by 1) genotype matrix
    // AF = (1^t K^(-1) 1)^(-1) * (1^t K^(-1) G)
    //    = (1^t U (lambda)^(-1) * U^t * 1)^(-1) *
    //      (1^t U (lambda)^(-1) * U^t * G)
    //    = ((U^t*1)^t * (lambda)^(-1) * (U^t*1))^(-1)
    //      ((U^t*1)^t * (lambda)^(-1) * (U^t * G))
    //    = (u1s * u1)^(-1) * (u1s * ug)
    const Eigen::MatrixXf& U = kinshipU.mat;
    Eigen::MatrixXf u1 = U.transpose().rowwise().sum();
    Eigen::ArrayXf u1s = u1.array() / (this->lambda.array()).abs().array();
    // This is usually the same, may cache the results to avoid future
    // repetitative calculations
    double denom = (u1s * u1.array()).sum();
    if (denom == 0.0) {
      return 0.0;
    }
    double numer = (u1s * ug.array()).sum();
    double beta = numer / denom;

    // here x is represented as 0, 1, 2, so beta(0, 0) is the mean genotype
    // multiply by 0.5 to get AF
    double af = 0.5 * beta;
    return af;
  }
  // NOTE: need to fit null model before calling this function
  // NOTE2: assume kinship matrices are unchanged
  double FastGetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                   Matrix& Xcol, int col) const {
    if (col >= Xcol.rows) return -1.;

    static bool initialized = false;
    static double denom;
    static Eigen::RowVectorXf alpha;

    const Eigen::MatrixXf& U = kinshipU.mat;
    if (!initialized) {
      Eigen::ArrayXf u1s;
      Eigen::MatrixXf u1 = U.transpose().rowwise().sum();  // n by 1
      u1s = u1.array() / (this->lambda.array()).abs().array();
      // denom = 1' * K^{-1} * 1 = u1' * lambda^{-1} * u1
      // here, we only use lambda (the kinship part), but not delta
      // aka. the genetic part (sigma_g), not the environment part (sigma_e)
      denom = (u1s * u1.array()).sum();
      // alpha = 1' * K^{-1}  = u1' * lambda^{-1} * u
      alpha = (u1.transpose() *
               this->lambda.array().abs().inverse().matrix().asDiagonal() *
               U.transpose())
                  .row(0);
      initialized = true;
      // fprintf(stderr, "initialized\n");
    }

    if (denom == 0.0) {
      return 0.0;
    }

    Eigen::MatrixXf g;
    G_to_Eigen(Xcol, &g);

    // double beta =  (alpha * g.array()).sum()  / denom;
    double beta = alpha.dot(g.col(col)) / denom;

    // here x is represented as 0, 1, 2, so beta(0, 0) is the mean genotype
    // multiply by 0.5 to get AF
    double af = 0.5 * beta;

    return af;
  }

  double GetPvalue() const { return this->pvalue; }
  double GetUStat() const { return this->Ustat; }
  double GetVStat() const { return this->Vstat; }
  double GetEffect() {
    if (this->Vstat != 0.0) return this->Ustat / this->Vstat;
    return 0.;
  };  // U/V
  double GetSE() {
    if (this->Vstat != 0.0) return 1 / sqrt(this->Vstat);
    return 0.;
  };  // 1/sqrt(V)
  double GetSigmaG2() { return this->sigma2; };
  double GetSigmaE2() { return this->sigma2 * this->delta; };
  double GetDelta() { return this->delta; };  // delta = sigma2_e / sigma2_g

  double GetNullLogLikelihood() const { return this->nullLikelihood; }
  double GetAltLogLikelihood() const { return this->altLikelihood; }
  void GetBeta(EigenMatrix* beta) const {
    EigenMatrix& b = *beta;
    b.mat = this->beta;
  }
  void GetNullCovEst(Vector* beta) {
    if (!beta) return;
    Eigen_Column_to_G(this->beta, 0, beta);
  }
  // Calculate (X' \Sigma^{-1} X)^{-1}
  // = (X'U (\simga_g2*\lambda + sigma_e2*I)U'X)^{-1}
  // = (X'U (\lambda + delta*I)U'X)^{-1} / \sigma_g2
  void GetNullCovB(Matrix* betaCov) {
    if (!betaCov) return;

    const int n = ux.cols();
    Eigen::MatrixXf m =
        (ux.transpose() * (lambda.array() + delta).matrix().asDiagonal() * ux)
            .ldlt()
            .solve(Eigen::MatrixXf::Identity(n, n));
    // fprintf(stderr, "m[0][0] = %g\n", m(0,0));
    Eigen_to_G(m, betaCov);
  }
  void CalculateFactors(const Eigen::MatrixXf& U) {
    const double sigmaE2 = sigma2 * delta;
    const double sigmaG2 = sigma2;

    sigmaK = sigmaG2 * lambda.array().square().sum() + sigmaE2 * lambda.sum();

    sigma1 =
        sigmaG2 * (U.colwise().sum().array().square().matrix() * lambda).sum() +
        sigmaE2 * U.rows();
  }
  double GetSigmaK() { return this->sigmaK; }
  double GetSigma1() { return this->sigma1; }
  // Z: covariate (U'Z stored in this->ux)
  // Y: phenotype/response (U'Y stored in this->uy)
  // Cov(ZY) = Z' Simga^{-1} Y
  //  = Z' (K*sigma_G2 + I*sigma_E2)^{-1} Y
  //  = Z' U (lambda*sigma_G2 + I*sigma_E2)^{-1} U' Y
  //  = (U'Z)' (lambda + I*delta)^{-1} (U' Y) / sigma_G2
  void GetCovZY(Matrix* zy) {
    Eigen::MatrixXf m;
    m = (ux.transpose() *
         (this->lambda.array() + delta).inverse().matrix().asDiagonal() * uy) /
        sigma2;
    Eigen_to_G(m, zy);
  }
  // similar to GetCovZY
  void GetCovZZ(Matrix* zz) {
    Eigen::MatrixXf m;
    m = (ux.transpose() *
         (this->lambda.array() + delta).inverse().matrix().asDiagonal() * ux) /
        sigma2;
    Eigen_to_G(m, zz);
  }
  void toEigen(const std::vector<double>& g, Eigen::MatrixXf* out) {
    const int n = g.size();
    (*out).resize(n, 1);
    for (int i = 0; i < n; ++i) {
      (*out)(i, 0) = g[i];
    }
  }
  void GetCovXX(const std::vector<double>& g1, const std::vector<double>& g2,
                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                double* out) {
    const int n = g1.size();
    Eigen::Map<const Eigen::MatrixXd> g1D(g1.data(), n, 1);
    Eigen::MatrixXf g1E = g1D.cast<float>();
    Eigen::Map<const Eigen::MatrixXd> g2D(g2.data(), n, 1);
    Eigen::MatrixXf g2E = g2D.cast<float>();

    *out =
        (g1E.array() * (this->lambda.array() + delta).inverse() * g2E.array())
            .sum() /
        this->sigma2;
  }
  void GetCovXX(FloatMatrixRef& g1, FloatMatrixRef& g2,
                const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                float* out) {
    // const int n = g1.size();
    // Eigen::Map<const Eigen::MatrixXd> g1D(g1.data(), n, 1);
    // Eigen::MatrixXf g1E = g1D.cast<float>();
    // Eigen::Map<const Eigen::MatrixXd> g2D(g2.data(), n, 1);
    // Eigen::MatrixXf g2E = g2D.cast<float>();

    REF_TO_EIGEN(g1, g1E);
    REF_TO_EIGEN(g2, g2E);
    *out =
        (g1E.array() * (this->lambda.array() + delta).inverse() * g2E.array())
            .sum() /
        this->sigma2;
  }
  void GetCovXZ(const std::vector<double>& g, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS, std::vector<double>* out) {
    // const Eigen::MatrixXf& U = kinshipU.mat;
    const int n = g.size();
    Eigen::Map<const Eigen::MatrixXd> gD(g.data(), n, 1);
    Eigen::MatrixXf gE = gD.cast<float>();

    // res: 1 by nCov matrix
    const int nCov = ux.cols();
    Eigen::MatrixXf res =
        gE.transpose() *
        (this->lambda.array() + delta).inverse().matrix().asDiagonal() * ux /
        sigma2;

    out->resize(nCov);
    for (int i = 0; i < nCov; ++i) {
      (*out)[i] = (double)res(0, i);
    }
  }
  void GetCovXZ(FloatMatrixRef& g, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS, FloatMatrixRef& out) {
    // // const Eigen::MatrixXf& U = kinshipU.mat;
    // const int n = g.size();
    // Eigen::Map<const Eigen::MatrixXd> gD(g.data(), n, 1);
    // Eigen::MatrixXf gE = gD.cast<float>();
    REF_TO_EIGEN(g, gE);
    REF_TO_EIGEN(out, gOut);

    // res: 1 by nCov matrix
    // const int nCov = ux.cols();
    assert(gOut.rows() == 1 && gOut.cols() == ux.cols());

    //    Eigen::MatrixXf res =
    gOut = gE.transpose() *
           (this->lambda.array() + delta).inverse().matrix().asDiagonal() * ux /
           sigma2;

    // out->resize(nCov);
    // for (int i = 0; i < nCov; ++i) {
    //   (*out)[i] = (double)res(0, i);
    // }
  }

  int TransformCentered(std::vector<double>* geno, const EigenMatrix& kinshipU,
                        const EigenMatrix& kinshipS) {
    // cast type double to float
    int n = geno->size();
    Eigen::Map<Eigen::MatrixXd> gD(geno->data(), n, 1);
    Eigen::MatrixXf g = gD.cast<float>();
    Eigen::RowVectorXf g_mean = g.colwise().mean();
    const Eigen::MatrixXf& U = kinshipU.mat;
    g = U.transpose() * (g.rowwise() - g_mean);

    // cast type back
    gD = g.cast<double>();
    return 0;
  }
  int TransformCentered(FloatMatrixRef& geno, const EigenMatrix& kinshipU,
                        const EigenMatrix& kinshipS) {
    // // cast type double to float
    // int n = geno->size();
    // Eigen::Map<Eigen::MatrixXd> gD(geno->data(), n, 1);
    // Eigen::MatrixXf g = gD.cast<float>();
    REF_TO_EIGEN(geno, g);
    Eigen::RowVectorXf g_mean = g.colwise().mean();
    const Eigen::MatrixXf& U = kinshipU.mat;
    g = U.transpose() * (g.rowwise() - g_mean);

    // // cast type back
    // gD = g.cast<double>();
    return 0;
  }
  int Transform(std::vector<double>* geno, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS) {
    // type conversion
    int n = geno->size();
    Eigen::Map<Eigen::MatrixXd> gD(geno->data(), n, 1);
    Eigen::MatrixXf g = gD.cast<float>();
    const Eigen::MatrixXf& U = kinshipU.mat;
    g = U.transpose() * g;

    // cast type back
    gD = g.cast<double>();
    return 0;
  }
  int Transform(FloatMatrixRef& geno, const EigenMatrix& kinshipU,
                const EigenMatrix& kinshipS) {
    // // type conversion
    // int n = geno->size();
    // Eigen::Map<Eigen::MatrixXd> gD(geno->data(), n, 1);
    // Eigen::MatrixXf g = gD.cast<float>();
    REF_TO_EIGEN(geno, g);
    const Eigen::MatrixXf& U = kinshipU.mat;
    g = U.transpose() * g;

    // // cast type back
    // gD = g.cast<double>();
    return 0;
  }
  int GetWeight(Vector* out) const {
    const int n = lambda.rows();
    out->Dimension(n);
    for (int i = 0; i < n; ++i) {
      (*out)[i] = sigma2 * (lambda(i) + delta);
    }
    // fprintf(stderr, "sigma2 = %g, lambda(0) = %g, lambda(99) = %g, delta =
    // %g\n", sigma2, lambda(0), lambda(99), delta);
    return 0;
  }
  void disableCenterGenotype() { this->needToCenterGentype = false; }

 private:
  // Eigen::MatrixXf S;
  float sigma2;  // sigma2_g
  float delta;   // delta =  sigma2_e / sigma2_g
  Eigen::MatrixXf beta;
  // temporary values
  Eigen::MatrixXf uy;  // U' * y
  Eigen::MatrixXf ux;  // U' * x
  Eigen::MatrixXf ug;  // U' * g
  Eigen::MatrixXf lambda;
  // for LRT test
  double nullLikelihood;
  double altLikelihood;
  // for score test
  Eigen::MatrixXf uResid;  // U' * (y - mean(y)) or U' * (y - X * beta)
  double Ustat;
  double Vstat;
  double stat;
  // for both tests
  double pvalue;

  FastLMM::Test test;
  FastLMM::Model model;
  double sumResidual2;  // sum (  (Uy - Ux *beta)^2/(lambda + delta) )

 private:
  Eigen::MatrixXf scaledK;
  double sigmaK;
  double sigma1;
  bool needToCenterGentype;

 private:
  static bool showDebug;
};

bool FastLMM::Impl::showDebug;

//////////////////////////////////////////////////
// FastLMM Interface
//////////////////////////////////////////////////
FastLMM::FastLMM(Test test, Model model) { this->impl = new Impl(test, model); }
FastLMM::~FastLMM() { delete this->impl; }

// @return 0 when success
int FastLMM::FitNullModel(Matrix& Xnull, Matrix& y, const EigenMatrix& kinshipU,
                          const EigenMatrix& kinshipS) {
  return this->impl->FitNullModel(Xnull, y, kinshipU, kinshipS);
}
int FastLMM::TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol,
                           const EigenMatrix& kinshipU,
                           const EigenMatrix& kinshipS) {
  return this->impl->TestCovariate(Xnull, y, Xcol, kinshipU, kinshipS);
}
int FastLMM::CalculateUandV(Matrix& Xnull, Matrix& Y, Matrix& Xcol,
                            const EigenMatrix& kinshipU,
                            const EigenMatrix& kinshipS, Matrix* uMat,
                            Matrix* vMat) {
  return this->impl->CalculateUandV(Xnull, Y, Xcol, kinshipU, kinshipS, uMat,
                                    vMat);
}
double FastLMM::GetAF(const EigenMatrix& kinshipU,
                      const EigenMatrix& kinshipS) {
  return this->impl->GetAF(kinshipU, kinshipS);
}
double FastLMM::GetAF(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                      Matrix& Xcol) {
  return this->impl->GetAF(kinshipU, kinshipS, Xcol);
}
double FastLMM::FastGetAF(const EigenMatrix& kinshipU,
                          const EigenMatrix& kinshipS, Matrix& Xcol) {
  return this->impl->FastGetAF(kinshipU, kinshipS, Xcol, 0);
}
double FastLMM::FastGetAF(const EigenMatrix& kinshipU,
                          const EigenMatrix& kinshipS, Matrix& Xcol, int col) {
  return this->impl->FastGetAF(kinshipU, kinshipS, Xcol, col);
}
double FastLMM::GetPvalue() { return this->impl->GetPvalue(); }
double FastLMM::GetUStat() const { return this->impl->GetUStat(); };
double FastLMM::GetVStat() const { return this->impl->GetVStat(); };
double FastLMM::GetEffect() const { return this->impl->GetEffect(); };  // U/V
double FastLMM::GetSE() const { return this->impl->GetSE(); };  // 1/sqrt(V)
double FastLMM::GetSigmaE2() const { return this->impl->GetSigmaE2(); };
double FastLMM::GetSigmaG2() const { return this->impl->GetSigmaG2(); };
double FastLMM::GetDelta() const {
  return this->impl->GetDelta();
};  // delta = sigma2_e / sigma2_g
double FastLMM::GetNullLogLikelihood() const {
  return this->impl->GetNullLogLikelihood();
};
double FastLMM::GetAltLogLikelihood() const {
  return this->impl->GetAltLogLikelihood();
};
void FastLMM::GetBeta(EigenMatrix* beta) const {
  return this->impl->GetBeta(beta);
}
void FastLMM::GetNullCovEst(Vector* beta) {
  return this->impl->GetNullCovEst(beta);
}
void FastLMM::GetNullCovB(Matrix* betaCov) {
  return this->impl->GetNullCovB(betaCov);
}
double FastLMM::GetSigmaK() { return this->impl->GetSigmaK(); }
double FastLMM::GetSigma1() { return this->impl->GetSigma1(); }
void FastLMM::GetCovZY(Matrix* zy) { return this->impl->GetCovZY(zy); }
void FastLMM::GetCovZZ(Matrix* zz) { return this->impl->GetCovZZ(zz); }
void FastLMM::GetCovXX(const std::vector<double>& g1,
                       const std::vector<double>& g2,
                       const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                       double* out) {
  return this->impl->GetCovXX(g1, g2, kinshipU, kinshipS, out);
}
void FastLMM::GetCovXX(FloatMatrixRef& g1, FloatMatrixRef& g2,
                       const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                       float* out) {
  return this->impl->GetCovXX(g1, g2, kinshipU, kinshipS, out);
}
void FastLMM::GetCovXZ(const std::vector<double>& g,
                       const EigenMatrix& kinshipU, const EigenMatrix& kinshipS,
                       std::vector<double>* out) {
  return this->impl->GetCovXZ(g, kinshipU, kinshipS, out);
}
void FastLMM::GetCovXZ(FloatMatrixRef& g, const EigenMatrix& kinshipU,
                       const EigenMatrix& kinshipS, FloatMatrixRef& out) {
  return this->impl->GetCovXZ(g, kinshipU, kinshipS, out);
}
int FastLMM::TransformCentered(std::vector<double>* x,
                               const EigenMatrix& kinshipU,
                               const EigenMatrix& kinshipS) {
  return this->impl->TransformCentered(x, kinshipU, kinshipS);
}
int FastLMM::TransformCentered(FloatMatrixRef& x, const EigenMatrix& kinshipU,
                               const EigenMatrix& kinshipS) {
  return this->impl->TransformCentered(x, kinshipU, kinshipS);
}
int FastLMM::Transform(std::vector<double>* x, const EigenMatrix& kinshipU,
                       const EigenMatrix& kinshipS) {
  return this->impl->Transform(x, kinshipU, kinshipS);
}
int FastLMM::Transform(FloatMatrixRef& x, const EigenMatrix& kinshipU,
                       const EigenMatrix& kinshipS) {
  return this->impl->Transform(x, kinshipU, kinshipS);
}
int FastLMM::GetWeight(Vector* out) const { return this->impl->GetWeight(out); }
void FastLMM::disableCenterGenotype() { this->impl->disableCenterGenotype(); }
// need to negaive the MLE to minize it
double goalFunction(double x, void* param) {
  FastLMM::Impl* p = (FastLMM::Impl*)param;
  p->getBetaSigma2(x);
  double ret = (-p->getLogLikelihood(x));
  return ret;
}
