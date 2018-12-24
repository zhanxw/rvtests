#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include "BoltPlinkLoader.h"

#include "regression/EigenMatrix.h"
#include "regression/EigenMatrixInterface.h"
#include "regression/MatrixOperation.h"
#include "third/eigen/Eigen/Dense"

#ifdef _OPENMP
// highly recommmended to include omp to speed up
#include <omp.h>
#else
static int omp_get_max_threads() { return 1; }
static int omp_get_thread_num() { return 0; }
#endif

const int BoltPlinkLoader::Mask[4] = {3, 3 << 2, 3 << 4, 3 << 6};
const int BoltPlinkLoader::Shift[4] = {0, 2, 4, 6};
// (HOM_REF, MISSING, HET, HOM_ALT) == (0, 1, 2, 3)
const float BoltPlinkLoader::Plink2Geno[4] = {0, -1, 1, 2};

BoltPlinkLoader::BoltPlinkLoader(CommonVariable& cv)
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
// BoltPlinkLoader(const std::string& fn) : pin_(fn) {
int BoltPlinkLoader::open(const std::string& fn) {
  pin_ = new PlinkInputFile(fn);
  M_ = pin_->getNumMarker();
  N_ = pin_->getNumSample();
  C_ = INT_MIN;
  BatchSize_ = 64;

  NumBatch_ = (M_ + BatchSize_) / BatchSize_;
  Nstride_ = (N_ + 3) / 4;  // bytes needed for one marker in PLINK
  M2_ = (M_ + BatchSize_) / BatchSize_ * BatchSize_;
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
BoltPlinkLoader::~BoltPlinkLoader() {
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
int BoltPlinkLoader::loadCovariate(const std::string& fn) {
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
int BoltPlinkLoader::extractCovariateBasis() {
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
int BoltPlinkLoader::preparePhenotype(const Matrix* phenotype,
                                      bool binaryMode) {
  y_ = Eigen::MatrixXf::Zero(N_ + C_, 1);
  if (phenotype) {
    for (int i = 0; i < (int)N_; ++i) {
      y_(i, 0) = (*phenotype)(i, 0);
    }
  } else {
    const std::vector<double>& pheno = pin_->getPheno();
    for (int i = 0; i < (int)N_; ++i) {
      y_(i, 0) = pheno[i];
    }
  }
  if (!binaryMode) {  // no need to center phenotype for binary trait
    float avg = y_.sum() / N_;
    y_.topLeftCorner(N_, 1).array() -= avg;
  }
  y_.bottomLeftCorner(C_, 1).noalias() =
      z_.transpose() * y_.topLeftCorner(N_, 1);

  return 0;
}
int BoltPlinkLoader::prepareGenotype() {
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
  std::vector<int> alleleCount(256, 0);   // num alt-allele counts
  std::vector<int> alleleCount2(256, 0);  // num of homAlt genotypes
  std::vector<int> missingCount(256, 0);  // num of missing genotypes
  for (int i = 0; i < 256; ++i) {
    for (int j = 0; j < 4; ++j) {
      int g = (i & Mask[j]) >> Shift[j];
      switch (g) {
        case PlinkInputFile::HET:
          alleleCount[i]++;
          // alleleCount2[i]++;
          break;
        case PlinkInputFile::HOM_ALT:
          alleleCount[i] += 2;
          // alleleCount2[i] += 4;
          alleleCount2[i]++;
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
    // const int numHomAlt = numAllele2;
    // const int numHet = numAllele - numAllele2 * 2;
    // const int numHomRef = N_ - numMissing - numHet - numHomAlt;
    const double af = 0.5 * numAllele / (N_ - numMissing);
    const double mean = af + af;
    // here we divide sqrt(2*p*q) as GCTA paper describes,
    // another normalization method is to divide by sqrt(sample variance)
    const double sd = sqrt(2.0 * af * (1.0 - af));
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

int BoltPlinkLoader::projectCovariate(Eigen::MatrixXf* mat) {
  assert(mat && (size_t)mat->rows() == N_ + C_);
  Eigen::MatrixXf& m = *mat;
  m.bottomRows(C_).noalias() = z_.transpose() * m.topRows(N_);
  return 0;
}

Eigen::MatrixXf BoltPlinkLoader::projectToCovariateSpace(
    const Eigen::MatrixXf& in) {
  assert(in.rows() == N_ + C_);
  Eigen::MatrixXf out = in.topRows(N_) - z_ * in.bottomRows(C_);
  return out;
}

Eigen::MatrixXf BoltPlinkLoader::predictedCovariateEffect() {
  return z_ * y_.bottomRows(C_);
}

// @param v is the any [N x 1] matrix, usually the predicted response of y
// @return y - v
Eigen::MatrixXf BoltPlinkLoader::getResidual(const Eigen::MatrixXf& v) {
  assert(v.rows() == N_);
  return y_.topRows(N_) - v;
}

void BoltPlinkLoader::buildTable(int m, float* table) {
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
int BoltPlinkLoader::loadSNPBatch(size_t batch, float* stage) {
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
int BoltPlinkLoader::loadSNPWithCovBatch(size_t batch, float* stage) {
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
int BoltPlinkLoader::loadSNPWithNegCovBatch(size_t batch, float* stage) {
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
int BoltPlinkLoader::loadRandomSNPWithCov(int nSnp, float* stage) {
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
    buildTable(indice[i], table);
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

    float* pCov = zg_.data() + indice[i] * zg_.rows();
    memcpy(p, pCov, sizeof(float) * C_);
  }

  // Eigen::Map<Eigen::MatrixXf> g(stage, (N_ + C_), nSnp);
  // g.bottomRows(C_).noalias() = z_.transpose() * g.topRows(N_);
  return 0;
}

const Eigen::MatrixXf& BoltPlinkLoader::getPhenotype() const { return y_; }
float* BoltPlinkLoader::getStage() const { return stage_; }
