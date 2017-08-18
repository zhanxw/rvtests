#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include "regression/BoltLMM.h"

#include "third/eigen/Eigen/Dense"
#include "third/gsl/include/gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q

#include "base/IO.h"
#include "base/TypeConversion.h"
#include "base/Utils.h"
#include "libVcf/PlinkInputFile.h"
#include "libsrc/MathMatrix.h"
#include "libsrc/Random.h"
#include "regression/EigenMatrix.h"
#include "regression/EigenMatrixInterface.h"
#include "regression/MatrixRef.h"

// 0: no debug info
// 1: some debug info
// 2: most debug info (timing)
// 3: most debug info (timing + intermediate file)
static int BOLTLMM_DEBUG = 0;

// #include <fstream>
// #include "base/Profiler.h"
#include "base/SimpleTimer.h"
#include "base/TimeUtil.h"

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

struct Float4 {
  float x[4];
  float operator[](int i) const { return x[i]; }
  float& operator[](int i) { return x[i]; }
};

// serveral modules share these settings
struct CommonVariable {
  CommonVariable() : random_(12345) {}
  size_t M_;
  size_t M2_;  // this is M round up to multiple of BatchSize_
  size_t N_;
  // size_t N2;       // this is N round up to multiple of 4
  size_t C_;  // >= 1 size_teger (size_tercept)
  size_t MCtrial_;
  size_t BatchSize_;  // batch size of marker are processed at a time
  size_t NumBatch_;   // number of batches
  Random random_;
};

class PlinkLoader {
 public:
  PlinkLoader(CommonVariable& cv)
      : pin_(NULL),
        M_(cv.M_),
        M2_(cv.M2_),
        N_(cv.N_),
        // N2_(cv.N2_),
        C_(cv.C_),
        BatchSize_(cv.BatchSize_),
        NumBatch_(cv.NumBatch_),
        random_(cv.random_),
        Nstride_(-1),
        genotype_(NULL),
        stage_(NULL),
        byte2genotype_(NULL) {}
  // PlinkLoader(const std::string& fn) : pin_(fn) {
  int open(const std::string& fn) {
    pin_ = new PlinkInputFile(fn);
    M_ = pin_->getNumMarker();
    N_ = pin_->getNumSample();
    C_ = INT_MIN;
    BatchSize_ = 64;

    NumBatch_ = (M_ + BatchSize_) / BatchSize_;
    Nstride_ = (N_ + 3) / 4;  // bytes needed for one marker in PLINK
    M2_ = (M_ + BatchSize_) / BatchSize_ * BatchSize_;
    // N2_ = (N_ + 4 ) / 4 * 4;
    genotype_ = NULL;  // defer memory allocation in prepareGenotype()
    snpLookupTable_.resize(M2_);
    stage_ = NULL;  // defer memory allocation until covariates are loaded
    if (omp_get_max_threads() > (int)BatchSize_) {
      fprintf(stderr,
              "Please specify OpenMP threads less than [ %d ], current value "
              "is [ %d ]\n",
              (int)BatchSize_, omp_get_max_threads());
      exit(1);
    }
    byte2genotype_ = new float[omp_get_max_threads() * 256 * 4 * sizeof(float)];

    return 0;
  }
  ~PlinkLoader() {
    if (genotype_) {
      delete[] genotype_;
      genotype_ = NULL;
    }
    if (stage_) {
      delete[] stage_;
      stage_ = NULL;
    }
    if (byte2genotype_) {
      delete[] byte2genotype_;
      byte2genotype_ = NULL;
    }
    delete pin_;
  }
  int loadCovariate(const std::string& fn) {
    FILE* fp = fopen(fn.c_str(), "r");
    if (fp == NULL) {
      // no .covar file => only intercept
      C_ = 1;
      z_ = Eigen::MatrixXf::Ones(N_, C_);
    } else {
      fclose(fp);
      C_ = 1;  // intercept
      const std::vector<std::string>& iid = pin_->getSampleName();
      LineReader lr(fn);
      std::vector<std::string> fd;
      int lineNo = 0;
      while (lr.readLineBySep(&fd, "\t ")) {
        ++lineNo;
        if (lineNo == 1) {
          if ((toupper(fd[0]) != "FID") || toupper(fd[1]) != "IID") {
            fprintf(stderr, "%s file does not have a proper header line!\n",
                    fn.c_str());
            return -1;
          }
          C_ += fd.size() - 2;  // first two columns are fid and iid
          z_.resize(N_, C_);
          continue;
        }
        if (fd[1] != iid[lineNo - 2]) {  // lineNo is 1-based and covar file
                                         // usually have header, so minus 2
          fprintf(
              stderr,
              "Mismatched order of PLINK FAM and covariate file on line %d!\n",
              lineNo);
          return (-1);
        }
        z_(lineNo - 2, 0) = 1.0;
        for (size_t i = 1; i != C_; ++i) {
          z_(lineNo - 2, i) = atof(fd[2 + i - 1]);
        }
      }
    }
    stage_ = new float[BatchSize_ * (N_ + C_)];
    return 0;
  }
  int extractCovariateBasis() {
    int numSingularValueKept = 1;
    if (N_ < 16) {
      Eigen::JacobiSVD<Eigen::MatrixXf> svd(z_, Eigen::ComputeThinU);
      float threshold = svd.singularValues()[0] * 1e-8;
      for (size_t i = 1; i != C_; ++i) {
        if (svd.singularValues()[i] > threshold) {
          numSingularValueKept++;
        }
      }
      z_ = svd.matrixU().leftCols(numSingularValueKept);
    } else {
      Eigen::BDCSVD<Eigen::MatrixXf> svd(z_, Eigen::ComputeThinU);
      float threshold = svd.singularValues()[0] * 1e-8;
      for (size_t i = 1; i != C_; ++i) {
        if (svd.singularValues()[i] > threshold) {
          numSingularValueKept++;
        }
      }
      z_ = svd.matrixU().leftCols(numSingularValueKept);
    }

    return 0;
  }
  int preparePhenotype(const Matrix* phenotype) {
    y_ = Eigen::MatrixXf::Zero(N_ + C_, 1);
    if (phenotype) {
      for (int i = 0; i < (int)N_; ++i) {
        y_(i, 0) = (*phenotype)[i][0];
      }
    } else {
      const std::vector<double>& pheno = pin_->getPheno();
      for (int i = 0; i < (int)N_; ++i) {
        y_(i, 0) = pheno[i];
      }
    }
    float avg = y_.sum() / N_;
    y_.array() -= avg;
    y_.bottomLeftCorner(C_, 1).noalias() =
        z_.transpose() * y_.topLeftCorner(N_, 1);
    return 0;
  }
  int prepareGenotype() {
    // load .bed file to the memory

    // convert to size_t is necessary
    // e.g. when Nstride_ * M2_ = 20000 * 123648 = 2.4 x 10^9
    // since the maximum 32bit integer is ~2 x 109,
    // there will be integer overflow, and thus crash the program with bad_alloc
    // NOTE: in gdb, this type of error can be caught using:
    //       b 'std::bad_alloc::bad_alloc()'
    // fprintf(stderr, "Allocate unsigned char [ %d * %d]", Nstride_, M2_);
    genotype_ = new unsigned char[Nstride_ * M2_];
    pin_->readBED(genotype_, M_ * Nstride_);

    // calculate maf, norms, build snpLookupTable
    std::vector<int> alleleCount(256, 0);
    std::vector<int> alleleCount2(256, 0);
    std::vector<int> missingCount(256, 0);
    for (int i = 0; i < 256; ++i) {
      for (int j = 0; j < 4; ++j) {
        int g = (i & Mask[j]) >> Shift[j];
        switch (g) {
          case PlinkInputFile::HET:
            alleleCount[i]++;
            alleleCount2[i]++;
            break;
          case PlinkInputFile::HOM_ALT:
            alleleCount[i] += 2;
            alleleCount2[i] += 4;
            break;
          case PlinkInputFile::MISSING:
            missingCount[i]++;
        }
      }
    }
    unsigned char* p = genotype_;
    for (size_t m = 0; m != M_; ++m) {
      assert(p == (genotype_ + m * Nstride_));
      int numAllele = 0;
      int numAllele2 = 0;
      int numMissing = 0;
      // NOTE: in PLINK when the number of samples are not multiple of 4
      // the remainder bits are 00 (homozygous REF)
      // since we count alternative alleles, we do not need to deal with the
      // remainder genotypes with special care.
      for (size_t n = 0; n != Nstride_; ++n) {
        numAllele += alleleCount[(*p)];
        numAllele2 += alleleCount2[(*p)];
        numMissing += missingCount[(*p)];
        ++p;
      }
      double mean = 1.0 * numAllele / (N_ - numMissing);
      double var = 1.0 *
                   (numAllele2 * (N_ - numMissing) - numAllele * numAllele) /
                   (N_ - numMissing) / (N_ - numMissing - 1);
      double sd = sqrt(var);
      if (sd > 0) {
        snpLookupTable_[m][PlinkInputFile::HOM_REF] = (0.0 - mean) / sd;
        snpLookupTable_[m][PlinkInputFile::HET] = (1.0 - mean) / sd;
        snpLookupTable_[m][PlinkInputFile::HOM_ALT] = (2.0 - mean) / sd;
        snpLookupTable_[m][PlinkInputFile::MISSING] = 0.0;
      } else {
        snpLookupTable_[m][PlinkInputFile::HOM_REF] = 0.0;
        snpLookupTable_[m][PlinkInputFile::HET] = 0.0;
        snpLookupTable_[m][PlinkInputFile::HOM_ALT] = 0.0;
        snpLookupTable_[m][PlinkInputFile::MISSING] = 0.0;
      }
    }
    for (size_t m = M_ + 1; m != M2_; ++m) {
      snpLookupTable_[m][PlinkInputFile::HOM_REF] = 0.0;
      snpLookupTable_[m][PlinkInputFile::HET] = 0.0;
      snpLookupTable_[m][PlinkInputFile::HOM_ALT] = 0.0;
      snpLookupTable_[m][PlinkInputFile::MISSING] = 0.0;
    }

    // prepare Z'X and gNorm2
    zg_.resize(C_, M2_);
    zg_.setZero();
    gNorm2_.resize(M2_);
    gNorm2_.setZero();

    // load SNP in batches
    for (size_t batch = 0; batch != NumBatch_; ++batch) {
      loadSNPBatch(batch, stage_);
      int lb = batch * BatchSize_;
      int ub = std::min(lb + BatchSize_, M_);

      // gBatch/g: [BatchSize_ x N_]
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_, ub - lb);
      zg_.block(0, lb, C_, ub - lb).noalias() = z_.transpose() * g;

      // calculate squared norm of (I-Z'Z)X
      gNorm2_.segment(lb, ub - lb) =
          g.colwise().squaredNorm() -
          zg_.block(0, lb, C_, ub - lb).colwise().squaredNorm();
    }
    return 0;
  }
  int projectCovariate(Eigen::MatrixXf* mat) {
    assert(mat && (size_t)mat->rows() == N_ + C_);
    Eigen::MatrixXf& m = *mat;
    m.bottomRows(C_).noalias() = z_.transpose() * m.topRows(N_);
    return 0;
  }

  void buildTable(int m, float* table) {
    const float* v = (float*)&snpLookupTable_[m];
    for (int i = 0; i < 256; ++i) {
      for (int j = 0; j < 4; ++j) {
        *(table++) = *(v + ((i & Mask[j]) >> Shift[j]));
      }
    }
  }
  // load data in batches
  // @param batch to gBatchSize_ [BatchSize_ * (N_)]
  //  @param, marker [batch*BatchSize_, min(
  // (batch+1)*BatchSize_, M_)]
  // will be extracted and stored in gBatchSize_
  int loadSNPBatch(size_t batch, float* stage) {
    // #ifdef DEBUG
    //     QuickTimer qt(__PRETTY_FUNCTION__);
    // #endif
    size_t lb = batch * BatchSize_;
    size_t ub = lb + BatchSize_;  // std::min(lb + BatchSize_, M_);
#pragma omp parallel for
    for (size_t i = lb; i < ub; ++i) {
      unsigned char* g = genotype_ + i * Nstride_;
      float* p = stage + (i - lb) * N_;
      // build lookup table
      float* table =
          byte2genotype_ + omp_get_thread_num() * 256 * 4 * sizeof(float);
      buildTable(i, table);
      // look up normalized genotypes
      const int strides = (N_ & ~0x3) >> 2;
      for (int j = 0; j < strides; ++j) {
        memcpy(p, table + *(g + j) * 4, sizeof(float) * 4);
        p += 4;
      }
      const int remainder = N_ & 0x3;
      if (remainder) {
        memcpy(p, table + *(g + Nstride_) * 4, sizeof(float) * remainder);
        p += remainder;
      }
      assert(p == stage + (i - lb) * N_ + (N_));
#if 0
      // naive method - slow
      for (int j = 0; j < N_; ++j) {
        const unsigned char gg = *(g + (j >> 2));
        const int offset = j & 0x3;
        *p = snpLookupTable_[i][(gg & Mask[offset]) >> Shift[offset]];
        p++;
      }
#endif
    }
    return 0;
  }
  // load batch @param batch to gBatchSize_ [BatchSize_ * (N_ + C_)]
  // for a given batch @param, marker [batch*BatchSize_, min(
  // (batch+1)*BatchSize_, M_)]
  // will be extracted and stored in gBatchSize_
  int loadSNPWithCovBatch(size_t batch, float* stage) {
    // #ifdef DEBUG
    //     QuickTimer qt(__PRETTY_FUNCTION__);
    // #endif
    size_t lb = batch * BatchSize_;
    size_t ub = lb + BatchSize_;  // std::min(lb + BatchSize_, M_);

#pragma omp parallel for
    for (size_t i = lb; i < ub; ++i) {
      unsigned char* g = genotype_ + i * Nstride_;
      float* p = stage + (i - lb) * (N_ + C_);
      // build lookup table
      float* table =
          byte2genotype_ + omp_get_thread_num() * 256 * 4 * sizeof(float);
      buildTable(i, table);
      // look up normalized genotypes
      const int strides = (N_ & ~0x3) >> 2;
      for (int j = 0; j < strides; ++j) {
        memcpy(p, table + *(g + j) * 4, sizeof(float) * 4);
        p += 4;
      }
      const int remainder = N_ & 0x3;
      if (remainder) {
        memcpy(p, table + *(g + Nstride_) * 4, sizeof(float) * remainder);
        p += remainder;
      }
#if 0
      for (int j = 0; j < N_; ++j) {
        const unsigned char gg = *(g + (j >> 2));
        const int offset = j & 0x3;
        *p = snpLookupTable_[i][(gg & Mask[offset]) >> Shift[offset]];
        p++;

      }
#endif
      assert(p == stage + (i - lb) * (N_ + C_) + (N_));

      float* pCov = zg_.data() + i * zg_.rows();
      memcpy(p, pCov, sizeof(float) * C_);

#if 0      
      for (int j = 0; j < C_; ++j) {
        *p = zg_(j, i);
        p++;
      }

#endif
    }
    return 0;
  }
  // load batch @param batch to gBatchSize_ [BatchSize_ * (N_ + C_)]
  // for a given batch @param, marker [batch*BatchSize_, min(
  // (batch+1)*BatchSize_, M_)]
  // will be extracted and stored in gBatchSize_
  int loadSNPWithNegCovBatch(size_t batch, float* stage) {
    // #ifdef DEBUG
    //     QuickTimer qt(__PRETTY_FUNCTION__);
    // #endif
    size_t lb = batch * BatchSize_;
    size_t ub = lb + BatchSize_;  // std::min(lb + BatchSize_, M_);
#pragma omp parallel for
    for (size_t i = lb; i < ub; ++i) {
      unsigned char* g = genotype_ + i * Nstride_;
      float* p = stage + (i - lb) * (N_ + C_);
      // build lookup table
      float* table =
          byte2genotype_ + omp_get_thread_num() * 256 * 4 * sizeof(float);
      buildTable(i, table);
      // look up normalized genotypes
      const int strides = (N_ & ~0x3) >> 2;
      for (int j = 0; j < strides; ++j) {
        memcpy(p, table + *(g + j) * 4, sizeof(float) * 4);
        p += 4;
      }
      const int remainder = N_ & 0x3;
      if (remainder) {
        memcpy(p, table + *(g + Nstride_) * 4, sizeof(float) * remainder);
        p += remainder;
      }
#if 0
      for (int j = 0; j < N_; ++j) {
        const unsigned char gg = *(g + (j >> 2));
        const int offset = j & 0x3;
        *p = snpLookupTable_[i][(gg & Mask[offset]) >> Shift[offset]];
        p++;
      }
#endif
      for (size_t j = 0; j != C_; ++j) {
        *p = -zg_(j, i);
        p++;
      }
      assert(p == stage + (i - lb) * (N_ + C_) + (N_ + C_));
    }
    return 0;
  }
  // load data in batches
  // @param nSnp into [nSnp * (N_+C_)]
  int loadRandomSNPWithCov(int nSnp, float* stage) {
    std::vector<size_t> indice(nSnp);
    for (int i = 0; i < nSnp; ++i) {
      indice[i] = (size_t)(random_.Next() * M_);
    }

// #ifdef DEBUG
//     QuickTimer qt(__PRETTY_FUNCTION__);
// #endif
#pragma omp parallel for
    for (int i = 0; i < nSnp; ++i) {
      unsigned char* g = genotype_ + indice[i] * Nstride_;
      float* p = stage + (i) * (N_ + C_);
      // build lookup table
      float* table =
          byte2genotype_ + omp_get_thread_num() * 256 * 4 * sizeof(float);
      buildTable(i, table);
      // look up normalized genotypes
      const int strides = (N_ & ~0x3) >> 2;
      for (int j = 0; j < strides; ++j) {
        memcpy(p, table + *(g + j) * 4, sizeof(float) * 4);
        p += 4;
      }
      const int remainder = N_ & 0x3;
      if (remainder) {
        memcpy(p, table + *(g + Nstride_) * 4, sizeof(float) * remainder);
        p += remainder;
      }
      assert(p == stage + (i) * (N_ + C_) + (N_));
#if 0      
      for (int j = 0; j < N_; ++j) {
        const unsigned char gg = *(g + (j >> 2));
        const int offset = j & 0x3;
        *p = snpLookupTable_[i][(gg & Mask[offset]) >> Shift[offset]];
        p++;
      }
#endif
    }

    Eigen::Map<Eigen::MatrixXf> g(stage, (N_ + C_), nSnp);
    g.bottomRows(C_).noalias() = z_.transpose() * g.topRows(N_);
    return 0;
  }

  const Eigen::MatrixXf& getPhenotype() const { return y_; }
  float* getStage() const { return stage_; }

 private:
  PlinkInputFile* pin_;

  // Common variables
  size_t& M_;
  size_t& M2_;  // this is M round up to multiple of BatchSize_
  size_t& N_;
  size_t& C_;
  size_t& BatchSize_;  // batch size of marker are processed at a time
  size_t& NumBatch_;   // number of batches
  Random& random_;

  size_t Nstride_;
  Eigen::MatrixXf z_;   // covariate [ N x C ] matrix
  Eigen::MatrixXf y_;   // phenotype [ N x C ] matrix
  Eigen::MatrixXf zg_;  // Z' * G , [ C x M ] matrix
  std::vector<Float4>
      snpLookupTable_;      // [M x 4] snpTable[i][j] store normalized values
  Eigen::VectorXf gNorm2_;  // vector norm of (g - Z Z' g)

  // the BED file content
  EIGEN_ALIGN16 unsigned char*
      genotype_;  // PLINK genotype matrix, SNP major, 2 bits/genotype

  // centered and scaled genotyped in a batch,
  // allocated to be a matrix of [(N+C) x BatchSize_ ]
  // used as [ N x BatchSize_ ] without covariate
  // or [(N+C) x BatchSize_ ] with covariates
  float* stage_;

  static const int Mask[4];
  static const int Shift[4];
  static const float Plink2Geno[4];
  // this is for fast loading genotypes
  // each OpenMP threads takes a 256 x 4 memory lot
  EIGEN_ALIGN16 float* byte2genotype_;
};  // class PlinkLoader

const int PlinkLoader::Mask[4] = {3, 3 << 2, 3 << 4, 3 << 6};
const int PlinkLoader::Shift[4] = {0, 2, 4, 6};
// (HOM_REF, MISSING, HET, HOM_ALT) == (0, 1, 2, 3)
const float PlinkLoader::Plink2Geno[4] = {0, -1, 1, 2};

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
  void init(PlinkLoader& pl) {
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
};

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
        NumBatch_(cv.NumBatch_) {}
  int FitNullModel(const std::string& prefix, Matrix* phenotype) {
    const char* boltLMMDebugEnv = std::getenv("BOLTLMM_DEBUG");
    if (boltLMMDebugEnv) {
      BoltLMM::BoltLMMImpl::showDebug = atoi(boltLMMDebugEnv);
    } else {
      BoltLMM::BoltLMMImpl::showDebug = 0;
    }

    // load phenotype, genotype
    std::string fn = prefix;
    if (pl.open(prefix)) {
      return -1;
    }
    // load covariates
    fn += ".covar";
    if (pl.loadCovariate(fn)) {
      fprintf(stderr, "Failed to load covariate file [ %s ]!\n", fn.c_str());
      return -1;
    }

    // normalize covariate
    pl.extractCovariateBasis();

    // regress out covariate from phenotype
    pl.preparePhenotype(phenotype);

    // project Z to G and
    // record mean and sd for each SNP
    pl.prepareGenotype();

    // get constants
    stage_ = pl.getStage();

    // calculate heritability
    EstimateHeritability();

    // calibrate a scaling factor
    EstimateInfStatCalibration(pl);

    // // calcualte alpha to speed up AF calculation in the future
    // Given (1) empirical kinship is not invertible; (2) N is large,
    // we just report the averaged allele frequencies.
    // CalculateAlpha();

    return 0;
  }

  // test @param Xcol
  int TestCovariate(Matrix& Xcol) {
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
    v_ = projNorm2(g_test_)(0) * H_inv_y_norm2_ * infStatCalibration_ / N_;

    if (v_ > 0.0) {
      effect_ = u_ / v_;
      pvalue_ = gsl_cdf_chisq_Q(u_ * u_ / v_, 1.0);
    } else {
      effect_ = 0.;
      pvalue_ = 1.0;
    }

    af_ = 0.5 * gg.sum() / gg.rows();

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

    double h2[7] = {0.};  // 7 is the limit of iterations used by BoltLMM
    double logDelta[7] = {0.};
    double f[7] = {0.};

#if 0
    // use MINQUE to get initial guess on h2
    // see: Zhou Xiang's bioarxiv paper
#ifdef DEBUG
    fprintf(stderr, "minque start - %s\n", currentTime().c_str());
#endif
    double S;
    if (g_.rows() < g_.cols()) {
      S = (g_ * g_.transpose()).squaredNorm() / M / M / (N - 1) /
          (N - 1) -
          1.0 / (N - 1);
    } else {
      S = (g_.transpose() * g_).squaredNorm() / M / M / (N - 1) /
          (N - 1) -
          1.0 / (N - 1);
    }
    double q = ((g_.transpose() * y_).squaredNorm() / M -
                (y_.squaredNorm())) /
        (N - 1) / (N - 1);
    double sigma2_g_est = q / S;
    double sigma2_e_g_est = y_.squaredNorm() / (N - 1);
    double minqueH2 = sigma2_g_est / sigma2_e_g_est;
#ifdef DEBUG
    fprintf(stderr, "minque S = %g\n", S);
    fprintf(stderr, "minque q = %g\n", q);
    fprintf(stderr, "minque sigma2_g_est = %g\n", sigma2_g_est);
    fprintf(stderr, "minque sigma2_e_est = %g\n",
            sigma2_e_g_est - sigma2_g_est);
    fprintf(stderr, "minque end - %s\n", currentTime().c_str());
#endif
#endif

    if (BoltLMM::BoltLMMImpl::showDebug >= 1) {
      fprintf(stderr, "bolt 1a) start - %s\n", currentTime().c_str());
    }
#if 0
    double initH2 = minqueH2;
    if (minqueH2 < 0 || minqueH2 > 1) {
      fprintf(stderr, "MINQUE h2 is out of range [ %g ]!\n", minqueH2);
      initH2 = 0.25;  // use BOLT-LMM init value
    }
#endif
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
      for (int b = 0; b < NumBatch_; ++b) {
        pl.loadSNPWithCovBatch(b, stage_);
        Eigen::Map<Eigen::MatrixXf> g_z(stage_, N_ + C_, BatchSize_);
        for (int i = 0; i < BatchSize_; ++i) {
          for (int j = 0; j < N_ + C_; ++j) {
            if (j) fputc('\t', fp);
            fprintf(fp, "%g", g_z(j, i));
          }
          fputc('\n', fp);
        }
      }
      fclose(fp);
      dumpToFile(w.y, "tmp.w.y");
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
  void solve(const Eigen::MatrixXf& y, const double delta,
             Eigen::MatrixXf* h_inv_y) {
    TIMER(__PRETTY_FUNCTION__);

    // 5e-4 is the default threshold used in BoltLMM paper
    // conjugateSolverBolt(g_, y, delta, 5e-4, h_inv_y);
    Eigen::MatrixXf& x = *h_inv_y;
    assert(x.rows() == y.rows() && y.cols() == y.cols());
    x.setZero();

    Eigen::MatrixXf r = y;
    Eigen::MatrixXf p = y;
    Eigen::VectorXf rsold;
    Eigen::VectorXf rsnew;
    Eigen::VectorXf ratio;  // rsnew / rsold
    Eigen::MatrixXf ap;     // [ (N+C) x (MCtrial+1) ]
    Eigen::VectorXf alpha;  // [ (MCtrial + 1) ]
    const int NUM_COL = y.cols();
    projNorm2(r, &rsold);

    const int MaxIter = 250;  // BOLT-LMM maximum iteration in conjugate solver
    const double Tol = 5e-4;  // BOLT-LMM tolerence
    const int maxIter = std::min((int)N_, MaxIter);
    for (int i = 0; i < maxIter; ++i) {
      computeHx(delta, p, &ap);
      alpha = rsold.array() / projDot(p, ap).array();
#ifdef DEBUG
      fprintf(stderr, "alpha:");
      for (int ii = 0; ii != NUM_COL; ++ii) {
        fprintf(stderr, "\t%f", alpha(ii));
      }
      fprintf(stderr, "\n");
#endif
      for (int ii = 0; ii != NUM_COL; ++ii) {
        if (!std::isfinite(alpha(ii)) || fabs(alpha(ii)) < Tol) {
          alpha(ii) = 0.0;
        }
      }
      x = x + p * alpha.asDiagonal();
      r = r - ap * alpha.asDiagonal();
      projNorm2(r, &rsnew);

      if ((rsnew.array() < Tol).all()) {
        break;
      }
      if ((rsnew - rsold).array().abs().maxCoeff() < Tol) {
        break;
      }
      ratio = rsnew.array() / rsold.array();
#ifdef DEBUG
      fprintf(stderr, "ratio:");
      for (int ii = 0; ii != NUM_COL; ++ii) {
        fprintf(stderr, "\t%f", ratio(ii));
      }
      fprintf(stderr, "\n");
      fprintf(stderr, "rsnew:");
      for (int ii = 0; ii != NUM_COL; ++ii) {
        fprintf(stderr, "\t%f", rsnew(ii));
      }
      fprintf(stderr, "\n");
      fprintf(stderr, "rsold:");
      for (int ii = 0; ii != NUM_COL; ++ii) {
        fprintf(stderr, "\t%f", rsold(ii));
      }
      fprintf(stderr, "\n");
#endif
      for (int ii = 0; ii != NUM_COL; ++ii) {
        if (rsnew(ii) < Tol || !std::isfinite(ratio(ii))) {
          ratio(ii) = 0.0;
        }
      }
      p = r + p * ratio.matrix().asDiagonal();

      rsold = rsnew;
      if (BoltLMM::BoltLMMImpl::showDebug >= 2) {
        fprintf(stderr, "i = %d\tdelta = %g", i, delta);
        for (int ii = 0; ii < rsnew.size(); ++ii) {
          fprintf(stderr, "\t%g", rsnew(ii));
        }
        fprintf(stderr, "\n");
        if ((rsnew.array() < -1.0).any()) {
          fprintf(stderr, "Norm2 should be always positive!\n");
        }
      }
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

  // calculate invser(K)*y, aka solve(K, y),
  // where K = G * G' / M
  // y: [ (N) x (1)]
  // NOTE: when K is GRM, K in singular and cannot be inverted
  void solveKinv(const Eigen::MatrixXf& y, Eigen::MatrixXf* k_inv_y) {
    TIMER(__PRETTY_FUNCTION__);
    // 5e-4 is the default threshold used in BoltLMM paper
    // conjugateSolverBolt(g_, y, delta, 5e-4, h_inv_y);
    Eigen::MatrixXf& x = *k_inv_y;
    assert(x.rows() == y.rows() && y.cols() == y.cols());
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

  // ret = H * y
  //   where H = X X' / M + delta * I
  // NOTE: X X' is [X; Z'X] * [X; -Z'X]'
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
    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithNegCovBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_ + C_, BatchSize_);
      X_y.block(b * BatchSize_, 0, BatchSize_, y.cols()).noalias() =
          g.transpose() * y;
    }

    for (size_t b = 0; b != NumBatch_; ++b) {
      pl.loadSNPWithCovBatch(b, stage_);
      Eigen::Map<Eigen::MatrixXf> g(stage_, N_ + C_, BatchSize_);
      (*ret).noalias() +=
          g * X_y.block(b * BatchSize_, 0, BatchSize_, y.cols());
    }
    const float invM = 1.0 / M_;
    (*ret) *= invM;
    (*ret).noalias() += delta * y;
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

  Eigen::VectorXf projDot(const Eigen::MatrixXf& v1,
                          const Eigen::MatrixXf& v2) {
    assert(C_ > 0);
    assert((size_t)v1.rows() == N_ + C_);
    assert((size_t)v2.rows() == N_ + C_);
    assert(v1.cols() == v2.cols());

    Eigen::VectorXf ret =
        (v1.topRows(N_).array() * v2.topRows(N_).array())
            .matrix()
            .colwise()
            .sum() -
        (v1.bottomRows(C_).array() * v2.bottomRows(C_).array())
            .matrix()
            .colwise()
            .sum();
    return ret;
  }
  void projDot(const Eigen::MatrixXf& v1, const Eigen::MatrixXf& v2,
               Eigen::VectorXf* ret) {
    assert(C_ > 0);
    assert((size_t)v1.rows() == N_ + C_);
    assert((size_t)v2.rows() == N_ + C_);
    assert(v1.cols() == v2.cols());

    (*ret).noalias() = (v1.topRows(N_).array() * v2.topRows(N_).array())
                           .matrix()
                           .colwise()
                           .sum() -
                       (v1.bottomRows(C_).array() * v2.bottomRows(C_).array())
                           .matrix()
                           .colwise()
                           .sum();
  }

  void projNorm2(const Eigen::MatrixXf& v, Eigen::VectorXf* ret) {
    projDot(v, v, ret);
  }
  Eigen::VectorXf projNorm2(const Eigen::MatrixXf& v) { return projDot(v, v); }

  // estimate scaling factor using some random SNPs
  int EstimateInfStatCalibration(PlinkLoader& pl) {
    TIMER(__PRETTY_FUNCTION__);

    int nSnp = BatchSize_;
    // BOLT-LMM uses 30
    nSnp = std::min(30, (int)M_);

    Eigen::MatrixXf g(N_ + C_, nSnp);
    pl.loadRandomSNPWithCov(nSnp, g.data());
    Eigen::MatrixXf V_inv_x(g.rows(), g.cols());
    solve(g, sigma2_g_est_, sigma2_e_est_, &V_inv_x);

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
      fprintf(stderr, "ratio = %f\n", x_V_inv_x.sum() / x_x.sum());
    }
    fprintf(stderr, "infStatCalibration_ = %f\n", infStatCalibration_);

    // calculate empirically X'HX / X'X
    xVx_xx_ratio_ = x_V_inv_x.sum() / x_x.sum();
    if (!std::isfinite(xVx_xx_ratio_)) {
      xVx_xx_ratio_ = 1.0;
    }
    fprintf(stderr, "\ni\tx_v_inv_x\tx_x\tratio\n");
    for (int i = 0; i < nSnp; ++i) {
      fprintf(stderr, "%d\t%f\t%f\t%f\n", i, x_V_inv_x(i), x_x(i),
              x_x(i) == 0 ? 0.0 : x_V_inv_x(i) / x_x(i));
    }
    fprintf(stderr, "ratio = %f\n", x_V_inv_x.sum() / x_x.sum());

    return 0;
  }

 private:
  CommonVariable cv;
  PlinkLoader pl;

  // Common variables
  size_t& M_;
  size_t& M2_;  // this is M round up to multiple of BatchSize_
  size_t& N_;
  size_t& C_;
  size_t& MCtrial_;
  size_t& BatchSize_;  // batch size of marker are processed at a time
  size_t& NumBatch_;   // number of batches

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
  double delta_;
  double h2_;
  Eigen::MatrixXf H_inv_y_;
  double H_inv_y_norm2_;
  double infStatCalibration_;
  double xVx_xx_ratio_;
  Eigen::MatrixXf g_test_;
  Eigen::MatrixXf alpha_;

 private:
  static int showDebug;
};

int BoltLMM::BoltLMMImpl::showDebug = 0;

//////////////////////////////////////////////////
// BoltLMM class
//////////////////////////////////////////////////
BoltLMM::BoltLMM() { impl_ = new BoltLMMImpl; }
BoltLMM::~BoltLMM() { delete impl_; }

int BoltLMM::FitNullModel(const std::string& prefix, Matrix* phenotype) {
  return impl_->FitNullModel(prefix, phenotype);
}
int BoltLMM::TestCovariate(Matrix& Xcol) { return impl_->TestCovariate(Xcol); }

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
