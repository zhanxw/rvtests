#ifndef _BOLTPLINKLOADER_H_
#define _BOLTPLINKLOADER_H_

#include <vector>

#include "base/MathMatrix.h"
#include "base/TypeConversion.h"
#include "base/Utils.h"
#include "libVcf/PlinkInputFile.h"
#include "libsrc/Random.h"
#include "regression/EigenMatrix.h"

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

class BoltPlinkLoader {
 public:
  BoltPlinkLoader(CommonVariable& cv);
  // BoltPlinkLoader(const std::string& fn) : pin_(fn) {
  int open(const std::string& fn);
  ~BoltPlinkLoader();
  int loadCovariate(const std::string& fn);
  int extractCovariateBasis();
  int preparePhenotype(const Matrix* phenotype, bool binaryMode);
  int prepareGenotype();

  int projectCovariate(Eigen::MatrixXf* mat);

  Eigen::MatrixXf projectToCovariateSpace(const Eigen::MatrixXf& in);

  Eigen::MatrixXf predictedCovariateEffect();

  // @param v is the any [N x 1] matrix, usually the predicted response of y
  // @return y - v
  Eigen::MatrixXf getResidual(const Eigen::MatrixXf& v);

  void buildTable(int m, float* table);

  // load data in batches
  // @param batch to gBatchSize_ [BatchSize_ * (N_)]
  //  @param, marker [batch*BatchSize_, min(
  // (batch+1)*BatchSize_, M_)]
  // will be extracted and stored in gBatchSize_
  int loadSNPBatch(size_t batch, float* stage);

  // load batch @param batch to gBatchSize_ [BatchSize_ * (N_ + C_)]
  // for a given batch @param, marker [batch*BatchSize_, min(
  // (batch+1)*BatchSize_, M_)]
  // will be extracted and stored in gBatchSize_
  int loadSNPWithCovBatch(size_t batch, float* stage);

  // load batch @param batch to gBatchSize_ [BatchSize_ * (N_ + C_)]
  // for a given batch @param, marker [batch*BatchSize_, min(
  // (batch+1)*BatchSize_, M_)]
  // will be extracted and stored in gBatchSize_
  int loadSNPWithNegCovBatch(size_t batch, float* stage);

  // load data in batches
  // @param nSnp into [nSnp * (N_+C_)]
  int loadRandomSNPWithCov(int nSnp, float* stage);

  const Eigen::MatrixXf& getPhenotype() const;
  float* getStage() const;

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
};  // class BoltPlinkLoader

#endif /* _BOLTPLINKLOADER_H_ */
