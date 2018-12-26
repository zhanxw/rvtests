#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include "regression/BoltLMM.h"
#include "regression/BoltPlinkLoader.h"

#include "third/cnpy/cnpy.h"
#include "third/eigen/Eigen/Dense"
#include "third/gsl/include/gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q

#include "base/IO.h"
#include "base/MathMatrix.h"
#include "base/TypeConversion.h"
#include "base/Utils.h"
#include "libVcf/PlinkInputFile.h"
#include "libsrc/Random.h"
#include "regression/EigenMatrix.h"
#include "regression/EigenMatrixInterface.h"
#include "regression/MatrixOperation.h"
#include "regression/MatrixRef.h"
#include "regression/SaddlePointApproximation.h"

// 0: no debug info
// 1: some debug info
// 2: most debug info (timing)
// 3: most debug info (timing + intermediate file)
static int BOLTLMM_DEBUG = 0;

// Helpful environment variables:
// BOLTLMM_SAVE_NULL_MODEL=filename => save null model
// BOLTLMM_LOAD_NULL_MODEL=filename => load null model, need to verify it loads
// correctly

// #define DEBUG

// #include <fstream>
// #include "base/Profiler.h"
#include "base/SimpleTimer.h"
#include "base/TimeUtil.h"

#if defined(__APPLE__) && !defined(_OPENMP)
static int omp_get_max_threads() { return 1; }
static int omp_get_thread_num() { return 0; }
#endif

class QuickTimer {
 public:
  QuickTimer(const std::string& msg) { start(msg); }
  ~QuickTimer() { stop(); }
  void start(const std::string& msg) {
    if (BOLTLMM_DEBUG >= 2) {
      msg_ = msg;
      fprintf(stderr, "%s: Enter %s \n", currentTime().c_str(), msg_.c_str());
      timer_.start();
    }
  }
  void stop() {
    if (BOLTLMM_DEBUG >= 2) {
      timer_.stop();
      fprintf(stderr, "%s: Exit %s - elapsed %g seconds \n",
              currentTime().c_str(), msg_.c_str(), timer_.getSeconds());
    }
  }

 private:
  AccurateTimer timer_;
  std::string msg_;
};

#define TIMER(msg) QuickTimer qt(msg);

class WorkingData {
 public:
  explicit WorkingData(CommonVariable& cv)
      : M_(cv.M_),
        M2_(cv.M2_),
        N_(cv.N_),
        C_(cv.C_),
        MCtrial_(cv.MCtrial_),
        BatchSize_(cv.BatchSize_),
        NumBatch_(cv.NumBatch_),
        random_(cv.random_),
        x_beta_rand(N_ + C_, MCtrial_),
        e_rand(N_ + C_, MCtrial_),

        H_inv_y(N_ + C_, MCtrial_ + 1),
        y(N_ + C_, MCtrial_ + 1) {
    // maybe use faster random number generator
    // see:
    // https://github.com/eddelbuettel/rcppziggurat/blob/master/man/ziggurat.Rd
  }
  void init(BoltPlinkLoader& pl) {
    // fill beta_rand with N(0, 1/M)
    // fill e_rand with N(0, 1)
    const double sqrtMInv = 1.0 / sqrt((double)M_);
    Eigen::MatrixXf beta_rand(M2_, MCtrial_);
    beta_rand.setZero();
    for (size_t i = 0; i != M_; ++i) {
      for (size_t j = 0; j != MCtrial_; ++j) {
        beta_rand(i, j) = random_.Normal() * sqrtMInv;
      }
    }
    // calculate X* beta
    x_beta_rand.setZero();
    float* stage_ = pl.getStage();
    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithCovBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_ + C_, BatchSize_);
      x_beta_rand.noalias() +=
          g * beta_rand.block(b * BatchSize_, 0, BatchSize_, MCtrial_);
    }

    for (size_t i = 0; i != N_; ++i) {
      for (size_t j = 0; j != MCtrial_; ++j) {
        e_rand(i, j) = random_.Normal();
      }
    }

    y.col(0) = pl.getPhenotype();
    pl.projectCovariate(&e_rand);

    assert((size_t)x_beta_rand.rows() == N_ + C_);
    assert((size_t)x_beta_rand.cols() == MCtrial_);
  }

 public:
  const size_t& M_;
  const size_t& M2_;  // this is M round up to multiple of BatchSize_
  const size_t& N_;
  const size_t& C_;
  const size_t& MCtrial_;
  const size_t& BatchSize_;  // batch size of marker are processed at a time
  const size_t& NumBatch_;   // number of batches
  Random& random_;

  Eigen::MatrixXf x_beta_rand;  // [(N+C) x (MCtrial)]
  Eigen::MatrixXf e_rand;       //  [(N+C) x (MCtrial)]

  // 1st column: y_data;
  // 2nd ... (MCtrial+1)columm: y_rand
  Eigen::MatrixXf H_inv_y;  // [(N+C) x (MCtrial+1)]
  Eigen::MatrixXf y;        // [(N+C) x (MCtrial+1)]
  // todo: add y_rand and y_obs for submatrices of y

  Eigen::MatrixXf beta_hat;  // [M2 x (MCtrial+1)]
  Eigen::MatrixXf e_hat;     //  [(N+C) x (MCtrial+1)]
};                           // class WorkingData

class BoltLMM::BoltLMMImpl {
 public:
  BoltLMMImpl()
      : pl(cv),
        M_(cv.M_),
        M2_(cv.M2_),
        N_(cv.N_),
        C_(cv.C_),
        MCtrial_(cv.MCtrial_),
        BatchSize_(cv.BatchSize_),
        NumBatch_(cv.NumBatch_),
        binaryMode(false),
        saddlePointCalculator(NULL),
        useSaddlePoint(false)  // default: turn off saddle point approximation
  {}
  virtual ~BoltLMMImpl() {
    if (saddlePointCalculator) {
      delete saddlePointCalculator;
    }
  }

  void enableBinaryMode() { this->binaryMode = true; }
  int FitNullModel(const std::string& prefix, const Matrix* phenotype) {
    const char* boltLMMDebugEnv = std::getenv("BOLTLMM_DEBUG");
    if (boltLMMDebugEnv) {
      BoltLMM::BoltLMMImpl::showDebug = atoi(boltLMMDebugEnv);
      fprintf(stderr, "BOLTLMM_DEBUG=%s\n", boltLMMDebugEnv);
    } else {
      BoltLMM::BoltLMMImpl::showDebug = 0;
    }

    // load phenotype, genotype
    std::string fn = prefix;
    if (pl.open(prefix)) {
      return -1;
    }
    // load covariates
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "Load covariate\n");
    }
    fn += ".covar";
    if (pl.loadCovariate(fn)) {
      fprintf(stderr, "Failed to load covariate file [ %s ]!\n", fn.c_str());
      return -1;
    }

    // normalize covariate
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "Normalize covariate\n");
    }
    pl.extractCovariateBasis();

    // regress out covariate from phenotype
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "Pre-process phenotypes\n");
    }
    pl.preparePhenotype(phenotype, binaryMode);

    // project Z to G and
    // record mean and sd for each SNP
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "Pre-process genotyeps\n");
    }
    pl.prepareGenotype();

    // get constants
    stage_ = pl.getStage();

    // load null model if necessary
    null_model_file_name_ = std::getenv("BOLTLMM_LOAD_NULL_MODEL")
                                ? std::getenv("BOLTLMM_LOAD_NULL_MODEL")
                                : "";
    if (null_model_file_name_.empty()) {
      fprintf(stderr, "Run BOLT Null Model\n");

      // calculate heritability
      EstimateHeritability();

      // calibrate a scaling factor
      EstimateInfStatCalibration(pl);
    } else {
      fprintf(stderr, "Load BOLT Null Model: %s\n",
              null_model_file_name_.c_str());
      cnpy::npz_t my_npz = cnpy::npz_load(null_model_file_name_);
      cnpy::NpyArray tmp = my_npz["H_inv_y"];
      Eigen::Map<Eigen::MatrixXf> tmp2(tmp.data<float>(), tmp.shape[0],
                                       tmp.shape[1]);
      H_inv_y_ = tmp2;
      H_inv_y_norm2_ = *my_npz["H_inv_y_norm2"].data<double>();
      infStatCalibration_ = *my_npz["infStatCalibration"].data<double>();
      xVx_xx_ratio_ = *my_npz["xVx_xx_ratio"].data<double>();
    }
    // // calcualte alpha to speed up AF calculation in the future
    // Given (1) empirical kinship is not invertible; (2) N is large,
    // we just report the averaged allele frequencies.
    // CalculateAlpha();

    // save null model
    null_model_file_name_ = std::getenv("BOLTLMM_SAVE_NULL_MODEL")
                                ? std::getenv("BOLTLMM_SAVE_NULL_MODEL")
                                : "";
    if (!null_model_file_name_.empty()) {
      fprintf(stderr, "Save BOLT Null Model: %s\n",
              null_model_file_name_.c_str());
      cnpy::npz_save(
          null_model_file_name_, "H_inv_y", H_inv_y_.data(),
          {(unsigned long)H_inv_y_.rows(), (unsigned long)H_inv_y_.cols()},
          "w");
      cnpy::npz_save(null_model_file_name_, "H_inv_y_norm2", &H_inv_y_norm2_,
                     {1}, "a");
      cnpy::npz_save(null_model_file_name_, "infStatCalibration",
                     &infStatCalibration_, {1}, "a");
      cnpy::npz_save(null_model_file_name_, "xVx_xx_ratio", &xVx_xx_ratio_, {1},
                     "a");
    }

    // saddle point approximation
    if (binaryMode && useSaddlePoint) {
      // printf("calculate saddle point\n");
      comptueBLUP(&mu_);
      y_ = pl.getPhenotype().topRows(N_);
      resid_ = pl.getResidual(mu_);
      saddlePointCalculator = new SaddlePointApproximation(y_, mu_, resid_);

      // prepare monte-carlo data
      saddlePointCalculator_mc_ = new SaddlePointApproximation*[10];
      resid_mc_.resize(10);
      y_mc_.resize(10);
      for (int i = 0; i < 10; ++i) {
        simulateMC(mu_, y_, &resid_mc_[i], &y_mc_[i]);
        saddlePointCalculator_mc_[i] =
            new SaddlePointApproximation(y_mc_[i], mu_, resid_mc_[i]);
      }

      fprintf(stderr, "y_mc = \n");
      for (int i = 0; i < 15; ++i) {
        for (int j = 0; j < 10; ++j) {
          fprintf(stderr, " %g", y_mc_[j](i, 0));
        }
        fprintf(stderr, "\n");
      }

      fprintf(stderr, "resid_mc = \n");
      for (int i = 0; i < 15; ++i) {
        for (int j = 0; j < 10; ++j) {
          fprintf(stderr, " %g", resid_mc_[j](i, 0));
        }
        fprintf(stderr, "\n");
      }
    }

    return 0;
  }  // end FitNullModel
  void simulateMC(const Eigen::MatrixXf& mu, const Eigen::MatrixXf& y,
                  Eigen::MatrixXf* resid_mc, Eigen::MatrixXf* y_mc) {
    (*resid_mc).resize(N_, 1);
    (*y_mc).resize(N_, 1);
    for (int i = 0; i < N_; ++i) {
      if (random_.Next() <= mu(i, 0)) {
        (*y_mc)(i, 0) = 1;
      } else {
        (*y_mc)(i, 0) = 0;
      }
    }
    (*resid_mc) = (*y_mc) - mu;
  }

  // test @param Xcol
  int TestCovariate(const Matrix& Xcol) {
    static Eigen::MatrixXf gg;
    G_to_Eigen(Xcol, &gg);
    if ((size_t)g_test_.rows() != N_ + C_) {
      g_test_.resize(N_ + C_, 1);
    }
    g_test_.topRows(N_) = gg;
    pl.projectCovariate(&g_test_);
    // u_ = (g_test_.transpose() * H_inv_y_)(0, 0);
    // v_ = (g_test_.squaredNorm() * H_inv_y_norm2_) * infStatCalibration_ /
    //      g_test_.rows();
    u_ = projDot(g_test_, H_inv_y_)(0, 0);
    // dumpToFile(H_inv_y_, "tmp.H_inv_y_");
    v_ = projNorm2(g_test_)(0) * H_inv_y_norm2_ * infStatCalibration_ / N_;

    if (v_ > 0.0) {
      effect_ = u_ / v_;
      pvalue_ = gsl_cdf_chisq_Q(u_ * u_ / v_, 1.0);
    } else {
      v_ = 0;  // reset to zero, this can occur due to float point arithmetic
      effect_ = 0.;
      pvalue_ = 1.0;
    }
    af_ = 0.5 * gg.sum() / gg.rows();

    // saddle point approximation
    if (binaryMode && useSaddlePoint) {
      if (pvalue_ < 1e-6) {
        fprintf(stderr, "Enable binary mode calibration via MonteCarlo\n");
        // calculate de-correlated x
        const Eigen::MatrixXf g_tilde = pl.projectToCovariateSpace(g_test_);
        Eigen::MatrixXf x_decorr;
        solve(g_tilde, delta_, &x_decorr);
        for (int i = 0; i < 10; ++i) {
          fprintf(stderr, "g_tilde(%d, 0) = %g\n", i, g_tilde(i, 0));
          fprintf(stderr, "x_decorr(%d, 0) = %g\n", i, x_decorr(i, 0));
        }

        // average p-values
        Eigen::MatrixXf p_value(10, 1);
        for (int i = 0; i < 10; ++i) {
          saddlePointCalculator_mc_[i]->calculatePvalue(x_decorr,
                                                        &p_value(i, 0));
        }

        float newPvalue = p_value.sum() / 10;
        if (true) {
          fprintf(stderr, "old_pvalue = %g, new_pvalue = %g, from [", pvalue_,
                  newPvalue);
          for (int i = 0; i < 10; ++i) {
            fprintf(stderr, "%g, ", p_value(i));
          }
          fprintf(stderr, "]\n");
        }
        if (std::isfinite(newPvalue)) {
          pvalue_ = newPvalue;
        }
      } else if (pvalue_ < 0.001) {
        fprintf(stderr, "Enable binary mode calibration\n");
        // binary mode correction
        const Eigen::MatrixXf g_tilde = pl.projectToCovariateSpace(g_test_);
        float newPvalue;
        if (!saddlePointCalculator->calculatePvalue(g_tilde, &newPvalue)) {
          fprintf(stderr, "old_pvalue = %g\t", pvalue_);
          fprintf(stderr, "new_pvalue = %g\n", newPvalue);
#ifdef DEBUG
          if (newPvalue * 100 < pvalue_) {
            fprintf(stderr, "Suspicious result!\n");
          }
#endif
          pvalue_ = newPvalue;
        }
      }
    }  // end if (binaryMode)
    return 0;
  }
#if 0
  void CalculateAlpha() {
    // // refer to GetAF() for formula details
    // const int N = g_.rows();
    Eigen::MatrixXf one = Eigen::MatrixXf::Ones(N_, 1);
    Eigen::MatrixXf k_inv_1(N_, 1);
    solveKinv(one, &k_inv_1); // => does not work well, GRM not invertable
    alpha_ = 0.5 / (one.transpose() * k_inv_1)(0, 0) * k_inv_1;
  }
#endif
  double GetAF() {
    // af = 0.5 * (1' * inv(K) * 1)^(-1) * (1' * inv(K) * g)
    //    = 0.5 * (1' * inv(K) * 1)^(-1) * (g' * inv(K) * 1)
    //    = g' * alpha in which
    //  alpha = 0.5 * (1' * inv(K) * 1)^(-1) * (inv(K) * 1)
    //  we can calculate inv(K) * 1 using conjugate gradient
    return af_;
  }
  double GetU() { return u_; }
  double GetV() { return v_; }
  double GetEffect() { return effect_; }
  double GetPvalue() { return pvalue_; }

  void GetCovXX(const std::vector<double>& g1, const std::vector<double>& g2,
                double* out) {
    assert(g1.size() == g2.size());
    assert(g1.size() == N_);

    Eigen::MatrixXf g(g1.size() + C_, 2);  //
    // copy in data
    for (size_t i = 0; i != N_; ++i) {
      g(i, 0) = g1[i];
      g(i, 1) = g2[i];
    }

    // project
    pl.projectCovariate(&g);

    // calculate
    (*out) = (double)projDot(g.col(0), g.col(1))(0);

    // scale
    (*out) *= xVx_xx_ratio_;
  }
  void GetCovXX(const FloatMatrixRef& g1, const FloatMatrixRef& g2,
                float* out) {
    assert(g1.nrow_ == g2.nrow_ && g1.ncol_ == g2.ncol_);
    assert(g1.nrow_ == N_);

    REF_TO_EIGEN(g1, g1E);
    REF_TO_EIGEN(g2, g2E);

    Eigen::MatrixXf g(N_ + C_, 2);  //
    // copy in data
    // for (size_t i = 0; i != N_; ++i) {
    //   g(i, 0) = g1[i];
    //   g(i, 1) = g2[i];
    // }
    g.col(0).head(N_) = g1E.col(0);
    g.col(1).head(N_) = g2E.col(0);

    // project
    pl.projectCovariate(&g);

    // calculate
    (*out) = projDot(g.col(0), g.col(1))(0);

    // scale
    (*out) *= xVx_xx_ratio_;
  }

 private:
  int EstimateHeritability() {
    TIMER(__PRETTY_FUNCTION__);
    MCtrial_ = std::max(std::min((int)(4e9 / N_ / N_), 15), 3);  // follow BOLT
    if (BOLTLMM_DEBUG >= 1) {
      fprintf(stderr, "MCtrial_ = %d\n", (int)MCtrial_);
    }
    WorkingData w(cv);
    w.init(pl);

    if (std::getenv("BOLTLMM_MINQUE")) {
      return EstimateHeritabilityMinque(&w);
    } else {
      return EstimateHeritabilityBolt(&w);
    }
  }
  int EstimateHeritabilityMinque(WorkingData* wd) {
    WorkingData& w = *wd;
    // w.x_beta_rand (N_+C_, MCtrial_);
    fprintf(stderr, "M_ = %d, N_ = %d, BS = %d\n", (int)M_, (int)N_,
            (int)BatchSize_);
    float* stage_ = pl.getStage();

    // dumpToFile(w.x_beta_rand, "tmp.x_beta_rand");
    // trace(K) = E(b' X X' b) = \sum( (X'b)_ij^2 ) , b is a column vector and
    // each b_i has 0 mean, unit variance
    float trace = (w.x_beta_rand.topRows(N_).squaredNorm() -
                   w.x_beta_rand.bottomRows(C_).squaredNorm()) /
                  MCtrial_;
    float trueTrace = 0;
    float trace2 = 0;  // trace(K^2)
    float gy = 0;
    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithCovBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_ + C_, BatchSize_);
      trueTrace += g.topRows(N_).squaredNorm() - g.bottomRows(C_).squaredNorm();
      trace2 += (g.topRows(N_).transpose() *  // N_ x BatchSize_
                 w.x_beta_rand.topRows(N_))
                    .squaredNorm() -
                (g.bottomRows(C_).transpose() *  // N_ x BatchSize_
                 w.x_beta_rand.bottomRows(C_))
                    .squaredNorm();
      ;
      gy += (w.y.col(0).head(N_).transpose() * g.topRows(N_)).squaredNorm() -
            (w.y.col(0).tail(C_).transpose() * g.bottomRows(C_)).squaredNorm();
    }
    trueTrace /= M_;
    trace2 = trace2 / MCtrial_ / M_;
    gy /= M_;
    float yNorm2 =
        w.y.col(0).head(N_).squaredNorm() - w.y.col(0).tail(C_).squaredNorm();

    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "empiricial trace(K) = %g\n", trace);
      fprintf(stderr, "randomized trace(K) = %g\n", trueTrace);
      fprintf(stderr, "theoretial trace(K) = %d\n", (int)(N_ - C_));
      fprintf(stderr, "randomized trace(K^2) = %g\n", trace2);
      fprintf(stderr, "y'Ky = %g \n", gy);
      fprintf(stderr, "y'y(yNorm2) = %g\n", yNorm2);
    }

#if 0
    // Xiang Zhou's MINQUE method may result in negative variance component estimation
    float q = (gy / M_ - yNorm2) / (N_ - 1) / (N_ - 1);
    float S = trace2 / (N_ - 1) / (N_ - 1) - 1.0 / (N_ - 1);
    double sigma2_g_est = q / S;
    double sigma2_e_g_est = (w.y.col(0).head(N_).squaredNorm() -
                             w.y.col(0).tail(C_).squaredNorm()) /
                            (N_ - 1);
    double minqueH2 = sigma2_g_est / sigma2_e_g_est;

    fprintf(stderr, "minque S = %g\n", S);
    fprintf(stderr, "minque q = %g\n", q);
    fprintf(stderr, "minque sigma2_g_est = %g\n", sigma2_g_est);
    fprintf(stderr, "minque sigma2_e_est = %g\n",
            sigma2_e_g_est - sigma2_g_est);
    fprintf(stderr, "minque h2 = %g\n", minqueH2);

    // these may be negative
    sigma2_g_est_ = sigma2_g_est;
    sigma2_e_est_ = sigma2_e_g_est - sigma2_g_est;
    h2_ = sigma2_g_est / sigma2_e_g_est;
#endif

    // use HE equation
    double sigma2_e_g_est = (w.y.col(0).head(N_).squaredNorm() -
                             w.y.col(0).tail(C_).squaredNorm()) /
                            (N_ - 1);
    // According to Haseman Elston regression:
    // \hat{sigma_g^2} = trace(yKy) / trace(K^2)
    sigma2_g_est_ = gy / trace2;
    sigma2_e_est_ = sigma2_e_g_est - sigma2_g_est_;
    h2_ = sigma2_g_est_ / sigma2_e_g_est;

    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "minque sigma2_g_est = %g\n", sigma2_g_est_);
      fprintf(stderr, "minque sigma2_e_est = %g\n", sigma2_e_est_);
      fprintf(stderr, "minque h2 = %g\n", h2_);
    }

    // compute Y_rand (Y_data is the first column)
    delta_ = sigma2_e_est_ / sigma2_g_est_;
    // computeY(w.x_beta_rand, w.e_rand, delta, &w.y);
    // solve H_inv_y, and beta
    solve(w.y.col(0), delta_, &w.H_inv_y);
    fprintf(stderr, "=> Finish solve(%d)\n", __LINE__);

    H_inv_y_ = w.H_inv_y.col(0) / sigma2_g_est_;
    H_inv_y_norm2_ = projNorm2(H_inv_y_)(0);

    return 0;
  }

  int EstimateHeritabilityBolt(WorkingData* wd) {
    WorkingData& w = *wd;
    double h2[7] = {0.};  // 7 is the limit of iterations used by BoltLMM
    double logDelta[7] = {0.};
    double f[7] = {0.};

    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "bolt 1a) start - %s\n", currentTime().c_str());
    }
    double initH2 = 0.25;  // use BOLT-LMM init value
    int i = 0;
    h2[i] = initH2;
    logDelta[i] = getLogDeltaFromH2(h2[i]);
    f[i] = evalREML(logDelta[i], w);
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "i = %d\tlogDelta = %f\tf = %f\tdelta = %f\th2 = %f\n", i,
              logDelta[i], f[i], exp(logDelta[i]), h2[i]);
    }

    i = 1;
    if (f[i - 1] < 0) {
      // h2[i] = 0.125;
      h2[i] = initH2 / 2;
    } else {
      // h2[i] = 0.5;
      h2[i] = std::min(initH2 * 2, 0.5 * initH2 + 0.5);
    }
    logDelta[i] = getLogDeltaFromH2(h2[i]);
    f[i] = evalREML(logDelta[i], w);
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "i = %d\tlogDelta = %f\tf = %f\tdelta = %f\th2 = %f\n", i,
              logDelta[i], f[i], exp(logDelta[i]), h2[i]);
    }

    for (i = 2; i < 7; ++i) {
      logDelta[i] = (logDelta[i - 2] * f[i - 1] - logDelta[i - 1] * f[i - 2]) /
                    (f[i - 1] - f[i - 2]);
      if (!std::isfinite(logDelta[i])) {
        --i;
        break;
      }
      // changed the upper threshold from 10 to 5
      // as it does not seem to be stable, and causes f(10) = nan
      if (logDelta[i] > 5) {
        logDelta[i] = 5;
      }
      if (logDelta[i] < -10) {
        logDelta[i] = -10;
      }
      h2[i] = 1. / (1. + exp(logDelta[i]));
      if (fabs(logDelta[i] - logDelta[i - 1]) < 0.01) {
        break;
      }
      f[i] = evalREML(logDelta[i], w);

      if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
        fprintf(stderr, "i = %d\tlogDelta = %f\tf = %f\tdelta = %f\th2 = %f\n",
                i, logDelta[i], f[i], exp(logDelta[i]), h2[i]);
      }
    }
    if (i == 7) {
      --i;
    }  // boudary case when i reaches its maximum
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      // print iterations
      fprintf(stderr, "\ni\tlogDelta\tf\tdelta\th2\n");
      for (int j = 0; j <= i; ++j) {
        double delta = exp(logDelta[j]);
        fprintf(stderr, "%d\t%f\t%f\t%f\t%f\n", j, logDelta[j], f[j], delta,
                h2[j]);
      }
      fprintf(stderr, "bolt 1a) end - %s\n", currentTime().c_str());
    }

    delta_ = exp(logDelta[i]);
    // sigma2_g_est_ = (y_.transpose() * w.H_inv_y_data)(0, 0) / (N_ - C_);
    sigma2_g_est_ = (projDot(w.y.col(0), w.H_inv_y.col(0)))(0) / (N_ - C_);
    sigma2_e_est_ = delta_ * sigma2_g_est_;
    assert(sigma2_g_est_ > 0);
    h2_ = h2[i];
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr,
              "i = %d, delta = %g, sigma2_g = %g, sigma2_e = %g, h2 = %g\n", i,
              delta_, sigma2_g_est_, sigma2_e_est_, h2_);
    }
    // need to multiply this
    // originally we calculate (K/M + delta)^(-1) * y
    // we now need (K/M * sigma2.g + delta * sigma2.e)^(-1) * y
    // so need to divide sigma2.g
    H_inv_y_ = w.H_inv_y.col(0) / sigma2_g_est_;
    H_inv_y_norm2_ = projNorm2(H_inv_y_)(0);
    return 0;
  }
  double evalREML(double logDelta, WorkingData& w) {
    TIMER(__PRETTY_FUNCTION__);

    double delta = exp(logDelta);

    if (BoltLMM::BoltLMMImpl::showDebug) {
      fprintf(stderr, "solve for data in evalREML()\n");
    }
    // compute Y_rand (Y_data is the first column)
    computeY(w.x_beta_rand, w.e_rand, delta, &w.y);
    // solve H_inv_y, and beta
    solve(w.y, delta, &w.H_inv_y);
    fprintf(stderr, "=> Finish solve(%d)\n", __LINE__);

    // w.beta_hat_data = g_.transpose() * w.H_inv_y_data / w.M;
    // w.e_hat_data = delta * w.H_inv_y_data;
    estimateBetaAndE(w.H_inv_y, delta, M_, &w.beta_hat, &w.e_hat);

    double r_rand[2] = {0.};  // numerator and denominator
    double r_data[2] = {0.};  // numerator and denominator

    r_data[0] = w.beta_hat.topLeftCorner(M_, 1).squaredNorm();
    r_data[1] = w.e_hat.topLeftCorner(N_, 1).squaredNorm() -
                w.e_hat.bottomLeftCorner(C_, 1).squaredNorm();
    r_rand[0] = w.beta_hat.topRightCorner(M_, MCtrial_).squaredNorm();
    r_rand[1] = w.e_hat.topRightCorner(N_, MCtrial_).squaredNorm() -
                w.e_hat.bottomRightCorner(C_, MCtrial_).squaredNorm();

    double f = log((r_data[0] / r_data[1]) / (r_rand[0] / r_rand[1]));
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "data = {%g, %g}, rand = {%g, %g}\n", r_data[0],
              r_data[1], r_rand[0], r_rand[1]);
    }
    if (BoltLMM::BoltLMMImpl::showDebug >= 3) {
      dumpToFile(w.y, "tmp.w.y");
      FILE* fp = fopen("tmp.g", "wt");
      for (size_t b = 0; b < NumBatch_; ++b) {
        pl.loadSNPWithCovBatch(b, stage_);
        Eigen::Map<Eigen::MatrixXf> g_z(stage_, N_ + C_, BatchSize_);
        for (size_t i = 0; i < BatchSize_; ++i) {
          for (size_t j = 0; j < N_ + C_; ++j) {
            if (j) fputc('\t', fp);
            fprintf(fp, "%g", g_z(j, i));
          }
          fputc('\n', fp);
        }
      }
      fclose(fp);
    }
    return f;

  }  // end evalREML
  // w.beta_hat_data = g_.transpose() * w.H_inv_y_data / w.M;
  // w.e_hat_data = delta * w.H_inv_y_data;

  // H_inv_y: [ (N+C) x (MCtrial+1)]
  // beta_hat: [ M x (MCtrial+1) ]
  // e_hat: [ (N+C) x (MCtrial+1)]
  void estimateBetaAndE(const Eigen::MatrixXf& H_inv_y, const double delta,
                        const int M, Eigen::MatrixXf* beta_hat,
                        Eigen::MatrixXf* e_hat) {
    TIMER(__PRETTY_FUNCTION__);

    // *beta_hat = g_.transpose() * H_inv_y / M;
    // Done batch load x
    (*beta_hat).resize(M2_, MCtrial_ + 1);
    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithCovBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g_z(stage_, N_ + C_, BatchSize_);
      (*beta_hat).block(b * BatchSize_, 0, BatchSize_, MCtrial_ + 1).noalias() =
          (g_z.topRows(N_).transpose() * H_inv_y.topRows(N_) -
           g_z.bottomRows(C_).transpose() * H_inv_y.bottomRows(C_)) /
          M_;
    }
    *e_hat = delta * H_inv_y;
  }

  // calculate invser(H)*y, aka solve(H, y),
  // where H = G * G' /M + delta * diag(N)
  // y: [ (N+C) x (MCtrial+1)] or [ (N+C) x (# of random SNPs) ]
  // delta: sigma2_e / sigma2_g
  void solve(const Eigen::MatrixXf& y, const double delta,
             Eigen::MatrixXf* h_inv_y) {
    TIMER(__PRETTY_FUNCTION__);
    fprintf(stderr, "=> Enter solve()\n");
    // 5e-4 is the default threshold used in BoltLMM paper
    // conjugateSolverBolt(g_, y, delta, 5e-4, h_inv_y);
    Eigen::MatrixXf& x = *h_inv_y;
    // x.resize(y.rows(), y.cols());
    // x.setZero();
    x = y.array() / delta;

    Eigen::MatrixXf r = y;  // [ (N+C) x (MCtrial + 1) ]
    Eigen::MatrixXf p = r;  // [ (N+C) x (MCtrial + 1) ]
    Eigen::RowVectorXf rsold;
    Eigen::RowVectorXf rsnew;
    Eigen::RowVectorXf ratio;  // rsnew / rsold
    Eigen::MatrixXf ap;        // [ (N+C) x (MCtrial+1) ]
    Eigen::RowVectorXf alpha;  // [ (MCtrial + 1) ]
    const int NUM_COL = y.cols();

    computeHx(delta, x, &ap);
    r = y - ap;
    p = r;
    projNorm2(r, &rsold);

    const int MaxIter = 250;  // BOLT-LMM maximum iteration in conjugate solver
    const double Tol = 5e-4;  // BOLT-LMM tolerence
    const int maxIter = std::min((int)N_, MaxIter);
    for (int i = 0; i < maxIter; ++i) {
      if (BOLTLMM_DEBUG >= 1) {
        fprintf(stderr, "i = %d, delta = %g\n", i, delta);
      }
      // char fn[50];
      computeHx(delta, p, &ap);
      // sprintf(fn, "tmp.cg.%d.p", i);
      // dumpToFile(p, fn);
      // sprintf(fn, "tmp.cg.%d.ap", i);
      // dumpToFile(ap, fn);
      // dumpToFile(rsold, "tmp.rsold");
      // dumpToFile(p, "tmp.p");
      // dumpToFile(ap, "tmp.ap");
      Eigen::RowVectorXf tmp = projDot(p, ap);
      if (BOLTLMM_DEBUG >= 1) {
        for (int ii = 0; ii < 5; ++ii) {
          fprintf(stderr, "p(%d) = %g\tap(%d) = %g\n", ii, p(ii), ii, ap(ii));
        }
      }
      alpha = rsold.array() / projDot(p, ap).array();
      // dumpToFile(alpha, "tmp.alpha");
      if (alpha.size() != NUM_COL) {
        fprintf(stderr, "%d != %d\n", (int)alpha.size(), NUM_COL);
        exit(1);
      }
#ifdef DEBUG
// fprintf(stderr, "alpha:");
// for (int ii = 0; ii != NUM_COL; ++ii) {
//   fprintf(stderr, "\t%.2g", alpha(ii));
// }
// fprintf(stderr, "\n");
#endif
      for (int ii = 0; ii != NUM_COL; ++ii) {
        // if (!std::isfinite(alpha(ii)) || fabs(alpha(ii)) < Tol) {
        if (!std::isfinite(alpha(ii))) {
          fprintf(stderr, "alpha(%d) = %g\n", ii, alpha(ii));
          alpha(ii) = 0.0;
        }
      }
      x = x + p * alpha.asDiagonal();
      r = r - ap * alpha.asDiagonal();
      projNorm2(r, &rsnew);

      if ((rsnew.array() < Tol).all()) {
        fprintf(stderr, "reach tolerence\n");
        break;
      }
      if ((rsnew - rsold).array().abs().maxCoeff() < Tol) {
        fprintf(stderr, "no improvement\n");
        break;
      }
      ratio = rsnew.array() / rsold.array();
#ifdef DEBUG
      fprintf(stderr, "\trsold\trsnew\talpha\tratio\n");
      for (int ii = 0; ii != NUM_COL; ++ii) {
        fprintf(stderr, "\t%.2g", rsold(ii));
        fprintf(stderr, "\t%.2g", rsnew(ii));
        fprintf(stderr, "\t%.2g", alpha(ii));
        fprintf(stderr, "\t%.2g", ratio(ii));
        fprintf(stderr, "\n");
      }
#endif
      for (int ii = 0; ii != NUM_COL; ++ii) {
        if (rsnew(ii) < Tol || !std::isfinite(ratio(ii))) {
          ratio(ii) = 0.0;
        }
      }
      p = r + p * ratio.asDiagonal();

      if (BoltLMM::BoltLMMImpl::showDebug >= 2) {
        // fprintf(stderr, "i = %d\tdelta = %g\tNorm2[", i, delta);
        // for (int ii = 0; ii < rsnew.size(); ++ii) {
        //   fprintf(stderr, "\t%.2g", rsnew(ii));
        // }
        // fprintf(stderr, "]\n");
        if ((rsnew.array() < -1.0).any()) {
          fprintf(stderr, "Norm2 should be always positive!\n");
        }
      }

      rsold = rsnew;
    }
  }

  // calculate invser(H)*y, aka solve(H, y),
  // where H = G * G' / M * sigma2_g + sigma2_e * diag(N)
  //         = sigma2_g * (G * G' / M * sigma2_g + delta * diag(N))
  //         = sigma2_g * HH
  // h_inv_y = inv(H) * y = inv(HH) * y / sigma2_g
  void solve(const Eigen::MatrixXf& y, const double sigma2_g,
             const double sigma2_e, Eigen::MatrixXf* h_inv_y) {
    assert(sigma2_g > 0.0);
    double delta = sigma2_e / sigma2_g;
    solve(y, delta, h_inv_y);  // h_inv_y is H^(-1)*y
    (*h_inv_y) /= sigma2_g;
  }
#if 0
  // calculate invser(K)*y, aka solve(K, y),
  // where K = G * G' / M
  // y: [ (N) x (1)]
  // NOTE: when K is GRM, K in singular and cannot be inverted, so do not use
  // the codes below.
  void solveKinv(const Eigen::MatrixXf& y, Eigen::MatrixXf* k_inv_y) {
    TIMER(__PRETTY_FUNCTION__);
    // 5e-4 is the default threshold used in BoltLMM paper
    // conjugateSolverBolt(g_, y, delta, 5e-4, h_inv_y);
    Eigen::MatrixXf& x = *k_inv_y;
    assert(x.rows() == y.rows() && x.cols() == y.cols());
    assert(y.cols() == 1);
    x.setZero();

    Eigen::MatrixXf r = y;
    Eigen::MatrixXf p = y;
    float rsold;
    float rsnew;
    Eigen::MatrixXf ap;  // [ (N) x (1) ]
    float alpha;         // [ (1) ]
    rsold = r.squaredNorm();

    const int MaxIter = 250;  // BOLT-LMM maximum iteration in conjugate solver
    const double Tol = 5e-4;  // BOLT-LMM tolerence
    const int maxIter = std::min((int)N_, MaxIter);
    for (int i = 0; i < maxIter; ++i) {
      computeKx(p, &ap);
      alpha = rsold / p.col(0).dot(ap.col(0));
      x = x + p * alpha;
      r = r - ap * alpha;
      rsnew = r.squaredNorm();

      if (rsnew < Tol) {
        break;
      }
      p = r + p * (rsnew / rsold);
      rsold = rsnew;

      if (BOLTLMM_DEBUG >= 1) {
        fprintf(stderr, "i = %d\talpha = %g\trsnew = %g\n", i, alpha, rsnew);
      }
    }
  }
#endif
  // ret = H * y
  //   where H = X X' / M + delta * I
  // In reality, H = X.plus * X.minus / M + delta * I
  //             y = [y; Z'y]
  //             ret = [XX'y - XX'ZZ'y; Z'(XX'y - XX'ZZ'y)]
  // NOTE: X X' is [X; Z'X] * [X; -Z'X]', denoted as X.plus and X.minus
  // respectively
  // y: [ (N+C) x (MCtrial + 1) ]
  void computeHx(double delta, const Eigen::MatrixXf& y, Eigen::MatrixXf* ret) {
    // #ifdef DEBUG
    //     QuickTimer qt(__PRETTY_FUNCTION__);
    // #endif
    // X: [X; Z'X] [ (N+C) x M ]
    //
    // X_y: [X' -X'Z] [ M by (MCtrial+1) ]
    ret->resize(N_ + C_, y.cols());
    ret->setZero();
    Eigen::MatrixXf X_y =
        Eigen::MatrixXf::Zero(M2_, y.cols());  // X' * y: [ M * (MCtrial + 1) ]

    // char fn[50];
    static int freq = 0;
    freq++;

    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithNegCovBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_ + C_, BatchSize_);
      X_y.block(b * BatchSize_, 0, BatchSize_, y.cols()).noalias() =
          g.transpose() * y;

      // if (freq <= 1) {
      //   sprintf(fn, "tmp.%d.snpNegCov.%d", freq, (int)b);
      //   dumpToFile(g, fn);
      // }
    }
    // sprintf(fn, "tmp.%d.X_y", freq);
    // dumpToFile(X_y, fn);

    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithCovBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_ + C_, BatchSize_);
      (*ret).noalias() +=
          g * X_y.block(b * BatchSize_, 0, BatchSize_, y.cols());

      // if (freq <= 1) {
      //   sprintf(fn, "tmp.%d.snpPosCov.%d", (int)b);
      //   dumpToFile(g, fn);
      // }
    }

    // sprintf(fn, "tmp.%d.ret", freq);
    // dumpToFile((*ret), fn);

    const float invM = 1.0 / M_;
    (*ret) *= invM;
    (*ret).noalias() += delta * y;

    // sprintf(fn, "tmp.Hx.%d.y", freq);
    // dumpToFile(y, fn);
    // sprintf(fn, "tmp.Hx.%d.ret", freq);
    // dumpToFile((*ret), fn);
  }

  // @param ret = BLUP(y) = X * X' / M * sigma2_g * H^(-1) * y + Z Z' y
  void comptueBLUP(Eigen::MatrixXf* ret) {
    // X: [X; Z'X] [ (N+C) x M ]
    //
    // X_y: [X' -X'Z] [ M by (MCtrial+1) ]
    ret->resize(N_, 1);
    ret->setZero();
    Eigen::MatrixXf X_y =
        Eigen::MatrixXf::Zero(N_, 1);  // X' * y: [ M * (MCtrial + 1) ]

    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPBatch(b, stage_);
      int lb = b * BatchSize_;
      int ub = std::min(lb + BatchSize_, M_);

      Eigen::Map<Eigen::MatrixXf> g(stage_, N_, ub - lb);
      X_y.block(lb, 0, ub - lb, 1).noalias() =
          g.transpose() * H_inv_y_.block(0, 0, N_, 1);
    }

    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithCovBatch(b, stage_);
      int lb = b * BatchSize_;
      int ub = std::min(lb + BatchSize_, M_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_, ub - lb);
      (*ret).noalias() += g * X_y.block(lb, 0, ub - lb, 1);
    }
    const float scale = sigma2_g_est_ / M_;
    (*ret) *= scale;
    (*ret).noalias() += pl.predictedCovariateEffect();
  }
  // ret = K * y
  //   where K = X X' / M
  // y: [ (N) x (1) ]
  void computeKx(const Eigen::MatrixXf& y, Eigen::MatrixXf* ret) {
    // #ifdef DEBUG
    //     QuickTimer qt(__PRETTY_FUNCTION__);
    // #endif
    // X: [X; Z'X] [ (N) x M ]
    //
    // X_y: [X'] [ M by (1) ]
    Eigen::MatrixXf X_y =
        Eigen::MatrixXf::Zero(M2_, y.cols());  // X' * y: [ M * (1) ]
    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_, BatchSize_);
      X_y.block(b * BatchSize_, 0, BatchSize_, y.cols()).noalias() =
          g.transpose() * y;
    }
    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_, BatchSize_);
      (*ret).noalias() +=
          g * X_y.block(b * BatchSize_, 0, BatchSize_, y.cols());
    }
    const float invM = 1.0 / M_;
    (*ret) *= invM;
    const float eps = 1e-3;
    (*ret).array() += eps;  // NOTE: has to add this to avoid convergence
                            // problem, as K can be ill-conditioned
  }
  // y_rand = g_ * beta_rand.col(tt) + sqrt(delta) * e_rand.col(tt);
  // beta: [ M2 x (MCtrial) ]
  // e: [ (N+C) x MCtrial ]
  // y: [ (N+C) x (MCtrial+1) ]
  void computeY(const Eigen::MatrixXf& x_beta, const Eigen::MatrixXf& e,
                const double delta, Eigen::MatrixXf* y) {
    // #ifdef DEBUG
    //     QuickTimer qt(__PRETTY_FUNCTION__);
    // #endif
    // *y = g_ * beta + sqrt(delta) * e;
    // Done: batch computation
    assert((size_t)x_beta.rows() == N_ + C_);
    assert((size_t)x_beta.cols() == MCtrial_);
    assert((size_t)e.rows() == N_ + C_);
    assert((size_t)e.cols() == MCtrial_);

    (*y).rightCols(MCtrial_).noalias() = x_beta + sqrt(delta) * e;
  }

  // delta = sigma2_e / sigma2_g (ref. eq 29)
  double getLogDeltaFromH2(const double h2) { return log((1.0 - h2) / h2); }

  Eigen::RowVectorXf projDot(const Eigen::MatrixXf& v1,
                             const Eigen::MatrixXf& v2) {
    assert(C_ > 0);
    assert((size_t)v1.rows() == N_ + C_);
    assert((size_t)v2.rows() == N_ + C_);
    assert(v1.cols() == v2.cols());

    // Eigen::VectorXf r1 = (v1.topRows(N_).array() * v2.topRows(N_).array())
    //                          .matrix()
    //                          .colwise()
    //                          .sum();
    // Eigen::VectorXf r2 = (v1.bottomRows(C_).array() *
    // v2.bottomRows(C_).array())
    //                          .matrix()
    //                          .colwise()
    //                          .sum();
    // dumpToFile(r1, "tmp.r1");
    // dumpToFile(r2, "tmp.r2");
    // Eigen::MatrixXf s1 = (v1.topRows(N_).array() *
    // v2.topRows(N_).array()).matrix();
    // Eigen::MatrixXf s2 = (v1.bottomRows(C_).array() *
    // v2.bottomRows(C_).array()).matrix();
    // dumpToFile(r1, "tmp.s1");
    // dumpToFile(r2, "tmp.s2");
    Eigen::RowVectorXf ret =
        (v1.topRows(N_).array() * v2.topRows(N_).array())
            .eval()
            .matrix()
            .colwise()
            .sum() -
        (v1.bottomRows(C_).array() * v2.bottomRows(C_).array())
            .eval()
            .matrix()
            .colwise()
            .sum();
    return ret;
  }
  void projDot(const Eigen::MatrixXf& v1, const Eigen::MatrixXf& v2,
               Eigen::RowVectorXf* ret) {
    assert(C_ > 0);
    assert((size_t)v1.rows() == N_ + C_);
    assert((size_t)v2.rows() == N_ + C_);
    assert(v1.cols() == v2.cols());

    (*ret).noalias() = (v1.topRows(N_).array() * v2.topRows(N_).array())
                           .eval()
                           .matrix()
                           .colwise()
                           .sum() -
                       (v1.bottomRows(C_).array() * v2.bottomRows(C_).array())
                           .eval()
                           .matrix()
                           .colwise()
                           .sum();
  }

  void projNorm2(const Eigen::MatrixXf& v, Eigen::RowVectorXf* ret) {
    assert(C_ > 0);
    assert((size_t)v.rows() == N_ + C_);

    *ret = (v.topRows(N_).cwiseAbs2().eval().colwise().sum() -
            v.bottomRows(C_).cwiseAbs2().eval().colwise().sum())
               .row(0);
    assert(v.cols() == ret->size());
  }

  Eigen::RowVectorXf projNorm2(const Eigen::MatrixXf& v) {
    assert(C_ > 0);
    assert((size_t)v.rows() == N_ + C_);

    Eigen::RowVectorXf ret =
        v.topRows(N_).cwiseAbs2().eval().colwise().sum() -
        v.bottomRows(C_).cwiseAbs2().eval().colwise().sum();
    return ret;
  }

  // estimate scaling factor using some random SNPs
  int EstimateInfStatCalibration(BoltPlinkLoader& pl) {
    TIMER(__PRETTY_FUNCTION__);

    int nSnp = BatchSize_;
    // BOLT-LMM uses 30
    nSnp = std::min(30, (int)M_);

    Eigen::MatrixXf g(N_ + C_, nSnp);
    pl.loadRandomSNPWithCov(nSnp, g.data());
    // dumpToFile(g, "tmp.g");
    Eigen::MatrixXf V_inv_x(g.rows(), g.cols());
    solve(g, sigma2_g_est_, sigma2_e_est_, &V_inv_x);
    // fprintf(stderr, "sigma2_g_est_ = %g, sigma2_e_est_ = %g\n",
    // sigma2_g_est_, sigma2_e_est_);
    // dumpToFile(V_inv_x, "tmp.V_inv_x");
    fprintf(stderr, "=> Finish solve(%d)\n", __LINE__);

    Eigen::VectorXf prospectiveStat(nSnp);
    Eigen::VectorXf uncalibratedRetrospectiveStat(nSnp);
    Eigen::VectorXf x_V_inv_y(nSnp);
    Eigen::VectorXf x_V_inv_x(nSnp);
    Eigen::VectorXf x_x;

    for (int i = 0; i < nSnp; ++i) {
      x_V_inv_y(i) = projDot(g.col(i), H_inv_y_)(0);
    }
    x_V_inv_x = projDot(g, V_inv_x);
    x_x = projNorm2(g);

    prospectiveStat = x_V_inv_y.array().square() / x_V_inv_x.array();
    uncalibratedRetrospectiveStat =
        (N_) * (x_V_inv_y.array().square()) / (x_x.array() * H_inv_y_norm2_);
    double r[2] = {0.};
    for (int i = 0; i < nSnp; ++i) {
      if (prospectiveStat(i) < 5.0) {
        r[0] += uncalibratedRetrospectiveStat(i);
        r[1] += prospectiveStat(i);
      }
    }
    if (r[1] != 0.0) {
      infStatCalibration_ = r[0] / r[1];
    } else {
      infStatCalibration_ = 1.0;
      fprintf(stderr, "Cannot estimate infStatCalibration_ [ %g / %g ]!\n",
              r[0], r[1]);
    }

    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr,
              "\ni\tprospectiveStat\tuncalibratedRetrospectiveStat\tratio\n");
      for (int i = 0; i < nSnp; ++i) {
        fprintf(stderr, "%d\t%f\t%f\t%f\n", i, prospectiveStat(i),
                uncalibratedRetrospectiveStat(i),
                prospectiveStat(i) / uncalibratedRetrospectiveStat(i));
      }
      fprintf(stderr, "infStatCalibration_ = %f\n", infStatCalibration_);
    }

    // calculate empirically X'HX / X'X
    xVx_xx_ratio_ = x_V_inv_x.sum() / x_x.sum();
    if (!std::isfinite(xVx_xx_ratio_)) {
      xVx_xx_ratio_ = 1.0;
    }
    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "\ni\tx_v_inv_x\tx_x\tratio\n");
      for (int i = 0; i < nSnp; ++i) {
        fprintf(stderr, "%d\t%f\t%f\t%f\n", i, x_V_inv_x(i), x_x(i),
                x_x(i) == 0 ? 0.0 : x_V_inv_x(i) / x_x(i));
      }
      fprintf(stderr, "ratio = %f\n", xVx_xx_ratio_);
    }

    return 0;
  }

 private:
  CommonVariable cv;
  BoltPlinkLoader pl;

  // Common variables
  const size_t& M_;
  const size_t& M2_;  // this is M round up to multiple of BatchSize_
  const size_t& N_;
  const size_t& C_;
  size_t& MCtrial_;
  const size_t& BatchSize_;  // batch size of marker are processed at a time
  const size_t& NumBatch_;   // number of batches

  double sigma2_g_est_;
  double sigma2_e_est_;

  // Eigen::MatrixXf
  //     g_;  // genotype matrix (Nsample by Mmarker) - too big to store
  EIGEN_ALIGN16 float*
      stage_;  // point to a memory for batch load genotype and covariate
  double af_;
  double u_;
  double v_;
  double effect_;
  double pvalue_;
  Random random_;
  double delta_;  //   sigma2_e / sigma2_g (ref. eq 29)
  double h2_;
  Eigen::MatrixXf H_inv_y_;  // H^(-1) * y , [(N+C) x 1] matrix
  double H_inv_y_norm2_;
  double infStatCalibration_;
  double xVx_xx_ratio_;
  Eigen::MatrixXf g_test_;
  // Eigen::MatrixXf alpha_;
  std::string null_model_file_name_;

  bool binaryMode;
  Eigen::MatrixXf mu_;     // BLUP of y = \sigma2_g * K * H_inv * y
  Eigen::MatrixXf y_;      // y [N x 1] matrix
  Eigen::MatrixXf resid_;  // resid = y - mu [N x 1] matrix

  SaddlePointApproximation* saddlePointCalculator;

  // Monte-carlo method
  std::vector<Eigen::MatrixXf> y_mc_;
  std::vector<Eigen::MatrixXf> resid_mc_;
  SaddlePointApproximation** saddlePointCalculator_mc_;
  bool useSaddlePoint;

 private:
  static int showDebug;
};  // end class BoltLMM::BoltLMMImpl

int BoltLMM::BoltLMMImpl::showDebug = 0;

//////////////////////////////////////////////////
// BoltLMM class
//////////////////////////////////////////////////
BoltLMM::BoltLMM() { impl_ = new BoltLMMImpl; }
BoltLMM::~BoltLMM() { delete impl_; }

int BoltLMM::FitNullModel(const std::string& prefix, const Matrix* phenotype) {
  return impl_->FitNullModel(prefix, phenotype);
}
int BoltLMM::TestCovariate(const Matrix& Xcol) {
  return impl_->TestCovariate(Xcol);
}

double BoltLMM::GetAF() { return impl_->GetAF(); }
double BoltLMM::GetU() { return impl_->GetU(); }
double BoltLMM::GetV() { return impl_->GetV(); }
double BoltLMM::GetEffect() { return impl_->GetEffect(); }
double BoltLMM::GetPvalue() { return impl_->GetPvalue(); }
void BoltLMM::GetCovXX(const std::vector<double>& g1,
                       const std::vector<double>& g2, double* out) {
  impl_->GetCovXX(g1, g2, out);
}
void BoltLMM::GetCovXX(const FloatMatrixRef& g1, const FloatMatrixRef& g2,
                       float* out) {
  impl_->GetCovXX(g1, g2, out);
}
void BoltLMM::enableBinaryMode() { impl_->enableBinaryMode(); }
//////////////////////////////////////////////////
// BoltLMM::BoltLMMImpl class
//////////////////////////////////////////////////

#if 0
// old code in BoltLMM::BoltLMMImpl
// not practical for large datasets due to exessive memory consumption
//
// P = I - Z * (Z' * Z)^(-1) * Z'
  void makeProjectionMatrix(Matrix& covar, Eigen::MatrixXf* ptr_proj_covar) {
    Eigen::MatrixXf& proj_covar = *ptr_proj_covar;

    if (covar.cols > 0) {
      Eigen::MatrixXf cov;
      G_to_Eigen(covar, &cov);

     // subtract intercept
      cov.rowwise() -= cov.colwise().mean();

      proj_covar = -cov * (cov.transpose() * cov).ldlt().solve(cov.transpose());
      // proj_covar = -cov * (cov.transpose() * cov).inverse() *
      // cov.transpose();
      proj_covar.diagonal().array() += 1;
      proj_covar.array() -= 1.0 / cov.rows();
    }
  }
  // P = I - Z * (Z' * Z)^(-1) * Z'
  void makeProjectionMatrix(int N, Eigen::MatrixXf* ptr_proj_covar) {
    Eigen::MatrixXf& proj_covar = *ptr_proj_covar;
    proj_covar.setConstant(N, N, -1.0 / N);
    proj_covar.diagonal().array() += 1;
  }
  void normalize(const Eigen::MatrixXf& genotype) {
    g_ = genotype;
    const int nrow = g_.rows();
    const int ncol = g_.cols();
    int numObs = 0;
    double sum = 0.;
    double sum2 = 0.;
    double avg;
    double sd_inv;
    for (int i = 0; i < ncol; ++i) {
      numObs = 0;
      sum = 0;
      sum2 = 0;
      for (int j = 0; j < nrow; ++j) {
        if (g_(j, i) < 0) {
          continue;
        }
        ++numObs;
        sum += g_(j, i);
        sum2 += g_(j, i) * g_(j, i);
      }
      if (numObs == 0) {
        avg = 0.0;
        sd_inv = 1.0;
      } else {
        avg = sum / numObs;
        sd_inv = 1.0 / sqrt(sum2 / numObs - avg * avg);
      }
      for (int j = 0; j < nrow; ++j) {
        if (g_(j, i) < 0) {
          g_(j, i) = 0;
        } else {
          g_(j, i) = (g_(j, i) - avg) * sd_inv;
        }
      }
    }
  }
#endif
