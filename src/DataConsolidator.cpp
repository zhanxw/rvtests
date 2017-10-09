#include "DataConsolidator.h"

#include <string>
#include <vector>

#include "base/Argument.h"
#include "base/CommonFunction.h"
#include "base/SimpleMatrix.h"
#include "libVcf/PlinkInputFile.h"
#include "libVcf/PlinkOutputFile.h"
#include "regression/EigenMatrix.h"
#include "regression/EigenMatrixInterface.h"
#include "src/LinearAlgebra.h"

DECLARE_BOOL_PARAMETER(boltPlinkNoCheck);

void convertToMinorAlleleCount(Matrix& in, Matrix* g) {
  Matrix& m = *g;
  m.Dimension(in.rows, in.cols);
  double s = 0;
  for (int j = 0; j < m.cols; ++j) {
    s = 0;
    for (int i = 0; i < m.rows; ++i) {
      s += in[i][j];
    }
    // maf = s / m.rows / 2
    // if maf > 0.5, s > m.row
    // otherwise, s <= m.row
    if (s <= m.rows) {
      for (int i = 0; i < m.rows; ++i) {
        m[i][j] = in[i][j];
      }
    } else {
      // flip to minor
      for (int i = 0; i < m.rows; ++i) {
        m[i][j] = 2 - in[i][j];
      }
    }
  }
}

void removeMissingMarker(Matrix* genotype) {
  Matrix& g = *genotype;
  // 1. find first nonmissing marker
  int missingCol = 0;
  for (missingCol = 0; missingCol < g.cols; ++missingCol) {
    if (hasMissingMarker(g, missingCol)) break;
  }
  if (missingCol == g.cols) return;  // no missing markers

  // 2. move non-missing sites to overwrite missing sites
  for (int col = missingCol + 1; col < g.cols; ++col) {
    if (hasMissingMarker(g, col)) {
      ++col;
    }

    for (int r = 0; r < g.rows; ++r) {
      g[r][missingCol] = g[r][col];
    }
    missingCol++;
  }
  g.Dimension(g.rows, missingCol);
}

bool isMonomorphicMarker(Matrix& genotype, int col) {
  if (col >= genotype.cols || col < 0) {
    logger->error("Invalid check of monomorhpic marker.");
    return false;
  }

  // find first non-missing genotype
  int nonMissingRow = genotype.rows;
  for (int i = 0; i < genotype.rows; ++i) {
    if (genotype[i][col] >= 0) {
      nonMissingRow = i;
      break;
    }
  }

  for (int r = nonMissingRow + 1; r < genotype.rows; ++r) {
    if (genotype[r][col] < 0)  // missing
      continue;
    if (genotype[r][col] != genotype[nonMissingRow][col]) return false;
  }
  return true;
}

void removeMonomorphicMarker(Matrix* genotype) {
  Matrix& g = *genotype;
  // 1. find first monomorhpic site
  int monoCol = 0;
  for (monoCol = 0; monoCol < g.cols; ++monoCol) {
    if (isMonomorphicMarker(g, monoCol)) {
      break;
    }
  }
  if (monoCol == g.cols) return;  // no monomorphic site

  // 2. move polymorphic sites to overwrite monomorphic sites
  for (int col = monoCol + 1; col < g.cols; ++col) {
    if (isMonomorphicMarker(g, col)) {
      continue;
    }

    // move g[, col] to g[, monoCol]
    for (int r = 0; r < g.rows; ++r) {
      g[r][monoCol] = g[r][col];
    }
    monoCol++;
  }
  g.Dimension(g.rows, monoCol);
}

DataConsolidator::DataConsolidator()
    : strategy(DataConsolidator::UNINITIALIZED),
      phenotypeUpdated(true),
      covariateUpdated(true),
      sex(NULL),
      formula(NULL),
      counter(NULL),
      parRegion(NULL) {}

DataConsolidator::~DataConsolidator() {}

/**
 * Impute missing genotype (<0) according to population frequency (p^2, 2pq,
 * q^2)
 */
void DataConsolidator::imputeGenotypeByFrequency(Matrix* genotype, Random* r) {
  Matrix& m = *genotype;
  for (int i = 0; i < m.cols; i++) {
    if ((*this->counter)[i].getNumMissing() == 0) {
      continue;
    }
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.rows; j++) {
      if (m[j][i] >= 0) {
        ac += m[j][i];
        an += 2;
      }
    }
    double p = an == 0 ? 0 : 1.0 * ac / an;
    double pRef = p * p;
    double pHet = pRef + 2.0 * p * (1.0 - p);
    for (int j = 0; j < m.rows; j++) {
      if (m[j][i] < 0) {
        double v = r->Next();
        if (v < pRef) {
          m[j][i] = 0;
        } else if (v < pHet) {
          m[j][i] = 1;
        } else {
          m[j][i] = 2;
        }
      }
    }
  }
}

/**
 * Impute missing genotype (<0) according to its mean genotype
 * @param genotype (people by marker matrix)
 */
void DataConsolidator::imputeGenotypeToMean(Matrix* genotype) {
  Matrix& m = *genotype;
  for (int i = 0; i < m.cols; i++) {
    if ((*this->counter)[i].getNumMissing() == 0) {
      continue;
    }
    int ac = 0;
    int an = 0;
    for (int j = 0; j < m.rows; j++) {
      if (m[j][i] >= 0) {
        ac += m[j][i];
        an += 2;
      }
    }
    double p;
    if (an == 0) {
      p = 0.0;
    } else {
      p = 1.0 * ac / an;
    }
    double g = 2.0 * p;
    for (int j = 0; j < m.rows; j++) {
      if (m[j][i] < 0) {
        m[j][i] = g;
      }
    }
    // fprintf(stderr, "impute to mean = %g, ac = %d, an = %d\n", g, ac, an);
  }
}

int DataConsolidator::preRegressionCheck(Matrix& pheno, Matrix& cov) {
  if (this->checkColinearity(cov)) {
    logger->warn(
        "The covariates is rank deficient and may suffer from colinearity!");
    return -1;
  }

  if (this->checkPredictor(pheno, cov)) {
    return -1;
  }
  return 0;
}

int DataConsolidator::checkColinearity(Matrix& cov) {
  // no covariate => no colienarity problem
  if (cov.cols == 0) {
    return 0;
  }

  const int r = matrixRank(cov);
  const int m = (cov.rows > cov.cols) ? cov.cols : cov.rows;
  if (m != r) {
    logger->warn(
        "The covariate matrix may be rank deficient (rank = [ %d ], dim = [ "
        "%d, %d ])!",
        r, cov.rows, cov.cols);
    return -1;
  }
  return 0;
}

int DataConsolidator::checkPredictor(Matrix& pheno, Matrix& cov) {
  Matrix& y = pheno;
  double r = 0.0;
  for (int i = 0; i < cov.cols; ++i) {
    if (corr(cov, i, y, 0, &r)) {
      logger->warn(
          "Failed to calculate correlation between covariate [ %s ] and "
          "phenotype",
          cov.GetColumnLabel(i));
      return -1;
    }
    if (fabs(r) > 0.999) {
      logger->warn(
          "Covariate [ %s ] has strong correlation [ r^2 = %g ] with the "
          "response!",
          cov.GetColumnLabel(i), r * r);
      return -1;
    }
  }
  return 0;
}

double DataConsolidator::getMarkerFrequency(int col) {
  return (*this->counter)[col].getAF();
}

void DataConsolidator::getMarkerFrequency(std::vector<double>* freq) {
  size_t n = (*this->counter).size();
  freq->resize(n);
  for (size_t i = 0; i != n; ++i) {
    (*freq)[i] = (*this->counter)[i].getAF();
  }
}

//////////////////////////////////////////////////
// codes related to kinship
int DataConsolidator::setKinshipSample(
    const std::vector<std::string>& samples) {
  if (samples.empty()) {
    logger->warn("There are no samples for the kinship");
    return -1;
  }
  this->kinship[KINSHIP_AUTO].setSample(samples);
  this->kinship[KINSHIP_X].setSample(samples);
  return 0;
}

int DataConsolidator::setKinshipFile(int kinshipType,
                                     const std::string& fileName) {
  if (kinshipType == KINSHIP_AUTO) {
    return this->kinship[kinshipType].setFile(fileName);
  }
  if (kinshipType == KINSHIP_X) {
    if (fileName.size()) {
      return this->kinship[KINSHIP_X].setFile(fileName);
    }
    // infer xHemi kinship file name
    std::string fn = this->kinship[KINSHIP_AUTO].getFileName();
    fn = fn.substr(0, fn.size() - 8);  // strip ".kinship"
    fn += ".xHemi.kinship";
    if (fileExists(fn)) {
      logger->info(
          "Kinship file [ %s ] detected and will be used for for X "
          "chromosome analysis",
          fn.c_str());
      return this->kinship[KINSHIP_X].setFile(fn);
    }
  }
  return -1;
}
int DataConsolidator::setKinshipEigenFile(int kinshipType,
                                          const std::string& fileName) {
  if (kinshipType == KINSHIP_AUTO) {
    return this->kinship[kinshipType].setEigenFile(fileName);
  }
  if (kinshipType == KINSHIP_X) {
    if (fileName.size()) {
      return this->kinship[KINSHIP_X].setEigenFile(fileName);
    }
    // infer xHemi kinship file name
    std::string fn = this->kinship[KINSHIP_AUTO].getEigenFileName();
    fn = fn.substr(0, fn.size() - 6);  // strip ".eigen"
    fn += ".xHemi.eigen";
    if (fileExists(fn)) {
      logger->info(
          "Kinship eigen-decomposition file [ %s ] detected and will be used "
          "for for X "
          "chromosome analysis",
          fn.c_str());
    }
    return this->kinship[KINSHIP_X].setEigenFile(fn);
  }
  return -1;
}
int DataConsolidator::loadKinship(int kinshipType) {
  return this->kinship[kinshipType].load();
}

int DataConsolidator::prepareBoltModel(
    const std::string& prefix, const std::vector<std::string>& sampleName,
    const SimpleMatrix& phenotype) {
  PlinkInputFile pin(prefix);
  if (!isUnique(sampleName)) {
    logger->error("%s:%d Unexpected duplicated samples!", __FILE__, __LINE__);
    exit(1);
  }
  for (size_t i = 0; i != sampleName.size(); ++i) {
    if (pin.getSampleIdx(sampleName[i]) < 0) {
      logger->warn(
          "%s:%d PLINK file [ %s ] does not include sample [ %s "
          "]!",
          __FILE__, __LINE__, (prefix + ".fam").c_str(), sampleName[i].c_str());
    }
  }

  const int M = pin.getNumMarker();
  const int N = pin.getNumIndv();
  bool needNewPlink = false;
  std::vector<int> sampleIdx;  // keep these samples
  std::vector<int> snpIdx;     // keep these SNPs
  if (FLAG_boltPlinkNoCheck) {
    // do nothing
    logger->info(
        "Assume binary PLINK file is ready for use (MAF, missingness)");
  } else {
    // calculate MAF
    std::vector<double> maf(M);
    pin.calculateMAF(&maf);

    // calculate missing rate
    std::vector<double> imiss(N);  // individual
    std::vector<double> lmiss(M);  // marker
    pin.calculateMissing(&imiss, &lmiss);
    // check missingness for samples
    std::set<int> badSampleIdx;
    for (int i = 0; i != N; ++i) {
      if (imiss[i] > 0.05) {
        logger->warn("Sample [ %s ] has high rate of missing genotype [ %g ]!",
                     pin.getIID()[i].c_str(), imiss[i]);
        badSampleIdx.insert(i);
        needNewPlink = true;
      }
    }
    if (badSampleIdx.size()) {
      logger->warn(
          "[ %d ] sample(s) have high missing rate, need to create new binary "
          "PLINK files",
          (int)badSampleIdx.size());
    }

    // choose SNPs to keep
    for (size_t i = 0; i != maf.size(); ++i) {
      if (maf[i] < 0.05 || lmiss[i] > 0.05) {
        needNewPlink = true;
        continue;
      } else {
        snpIdx.push_back(i);
      }
    }
    if ((int)snpIdx.size() != M) {
      logger->warn(
          "[ %d ] markers have high missing rate or low MAF, need to create "
          "new "
          "binary PLINK files",
          M - (int)snpIdx.size());
    }

    // build a sample index, such that plink.fam[index] is in the same order as
    // @param sampleName
    for (size_t i = 0; i != sampleName.size(); ++i) {
      if (badSampleIdx.count(i)) {
        continue;
      }
      sampleIdx.push_back(pin.getSampleIdx(sampleName[i]));
    }
    // check order
    for (size_t i = 0; i != sampleName.size(); ++i) {
      if (pin.getSampleIdx(sampleName[i]) != (int)i) {
        needNewPlink = true;
        logger->warn(
            "To adjust phenotype order, need to create new binary PLINK files");
        break;
      }
    }
    if ((int)sampleName.size() != N) {
      logger->warn(
          "Existing binary PLINK file has more samples [ %d ], need to create "
          "new binary PLINK files",
          N - (int)sampleName.size());
      needNewPlink = true;
    }
  }  // end    if (FLAG_boltPlinkNoCheck) {
  if (needNewPlink) {
    // write a new set of PLINK file
    this->boltPrefix = prefix + ".out";
    logger->info("Need to create new binary PLINK files with prefix [ %s ].",
                 this->boltPrefix.c_str());
    PlinkOutputFile pout(prefix + ".out");
    pout.extractFAMWithPhenotype(pin, sampleIdx, phenotype);
    pout.extractBIM(pin, snpIdx);
    pout.extractBED(pin, sampleIdx, snpIdx);
  } else {
    this->boltPrefix = prefix;
    logger->info("Use binary PLINK file [ %s ] in BoltLMM model",
                 this->boltPrefix.c_str());
  }

  // write covariate file if there are covariates
  if (covariate.cols > 0) {
    std::string covarFn = this->boltPrefix + ".covar";
    logger->info("Create covariate file [ %s ]", covarFn.c_str());
    FILE* fpCov = fopen((covarFn).c_str(), "wt");
    fprintf(fpCov, "FID\tIID");
    for (int j = 0; j < covariate.cols; ++j) {
      fprintf(fpCov, "\t%s", (char*)covariate.GetColumnLabel(j));
    }
    fprintf(fpCov, "\n");
    for (int i = 0; i < covariate.rows; ++i) {
      fprintf(fpCov, "%s\t%s", sampleName[i].c_str(), sampleName[i].c_str());
      for (int j = 0; j < covariate.cols; ++j) {
        fprintf(fpCov, "\t%g", covariate[i][j]);
      }
      fputs("\n", fpCov);
    }
    fclose(fpCov);
  }

  return 0;
}

#if 0
int DataConsolidator::loadGenotype(const std::string& prefix) {
  std::string fn = prefix;
  fn += ".dim";

  int nrow = 0;
  int ncol = 0;
  FILE* fp = fopen(fn.c_str(), "rt");
  fscanf(fp, "%d", &nrow);
  fscanf(fp, "%d", &ncol);
  fclose(fp);

  fullGenotype_ = new EigenMatrix;
  Eigen::MatrixXf& geno = fullGenotype_->mat;
  geno.resize(nrow, ncol);

  fn = prefix;
  fn += ".data";
  fp = fopen(fn.c_str(), "rb");
  double* buff = new double[nrow];

  for (int i = 0; i < ncol; i++) {
    // process column by column
    fread(buff, sizeof(double), nrow, fp);
    for (int j = 0; j < nrow; ++j) {
      geno(j, i) = buff[j];
    }
  }
  delete[] buff;
  fclose(fp);

  return 0;
}

int DataConsolidator::loadNormalizedGenotype(const std::string& prefix) {
  std::string fn = prefix;
  fn += ".dim";

  int nrow = 0;
  int ncol = 0;
  FILE* fp = fopen(fn.c_str(), "rt");
  fscanf(fp, "%d", &nrow);
  fscanf(fp, "%d", &ncol);
  fclose(fp);

  fullGenotype_ = new EigenMatrix;
  Eigen::MatrixXf& geno = fullGenotype_->mat;
  geno.resize(nrow, ncol);

  fn = prefix;
  fn += ".data";
  fp = fopen(fn.c_str(), "rb");
  double* buff = new double[nrow];
  double sum = 0;
  double sum2 = 0;
  int nonMiss = 0;

  double avg, sdInv;
  for (int i = 0; i < ncol; i++) {
    // process column by column
    fread(buff, sizeof(double), nrow, fp);
    sum = 0;
    sum2 = 0;
    nonMiss = 0;
    for (int j = 0; j < nrow; ++j) {
      if (buff[j] < 0) {
        continue;
      }
      geno(j, i) = buff[j];
      sum += buff[j];
      sum2 += buff[j] * buff[j];
      nonMiss++;
    }

    if (nonMiss == 0) {
      avg = 0.0;
      sdInv = 1.0;
    } else {
      avg = sum / nonMiss;
      sdInv = 1.0 / sqrt(sum2 / nonMiss - avg * avg);
    }

    for (int j = 0; j < nrow; ++j) {
      if (buff[j] < 0) {
        geno(j, i) = 0.;
      } else {
        geno(j, i) = (geno(j, i) - avg) * sdInv;
      }
    }
  }
  delete[] buff;
  fclose(fp);

  return 0;
}

EigenMatrix* DataConsolidator::getFullGenotype() { return fullGenotype_; }
#endif
