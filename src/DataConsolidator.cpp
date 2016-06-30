#include "DataConsolidator.h"

#include <string>
#include <vector>

#include "base/CommonFunction.h"
#include "regression/EigenMatrixInterface.h"

#include "LinearAlgebra.h"

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
};

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
};

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
    logger->warn("The covariate matrix may be rank deficient (%d < %d)!", r, m);
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
