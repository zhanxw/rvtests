#include "DataConsolidator.h"

#include <string>
#include <vector>

#include "base/CommonFunction.h"
#include "regression/EigenMatrixInterface.h"

#include "LinearAlgebra.h"
#include "snp_hwe.c"

void convertToMinorAlleleCount(Matrix& in, Matrix* g) {
  Matrix& m = *g;
  m.Dimension(in.rows, in.cols);
  double s = 0;
  for (int j = 0; j < m.cols; ++j) {
    s = 0;
    for (int i = 0; i < m.rows; ++i) {
      s += in[i][j];
    }
    if (2.0 * s < m.rows) {
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

double GenotypeCounter::getHWE() const {
  double hweP = 0.0;
  if (nHomRef + nHet + nHomAlt == 0 ||
      (nHet < 0 || nHomRef < 0 || nHomAlt < 0)) {
    hweP = 0.0;
  } else {
    hweP = SNPHWE(nHet, nHomRef, nHomAlt);
  }
  return hweP;
}

DataConsolidator::DataConsolidator()
    : strategy(DataConsolidator::UNINITIALIZED),
      phenotypeUpdated(true),
      covariateUpdated(true),
      sex(NULL),
      parRegion(NULL) {}

DataConsolidator::~DataConsolidator() {}

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
  int r = matrixRank(cov);
  int m = (cov.rows > cov.cols) ? cov.cols : cov.rows;
  // fprintf(stderr, "rank of cov is %d, and min(r,c) = %d\n", r, m);
  if (m != r) {
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
