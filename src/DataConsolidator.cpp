#include "DataConsolidator.h"

#include "third/eigen/Eigen/Core"
#include "third/eigen/Eigen/Eigenvalues"

#include "CommonFunction.h"
#include "regression/EigenMatrix.h"
#include "IO.h"

#undef DEBUG
//#define DEBUG
#ifdef DEBUG
#include <iostream>
#include <fstream>
#endif

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
  for (missingCol = 0; missingCol < g.cols; ++ missingCol) {
    if (hasMissingMarker(g, missingCol))
      break;
  }
  if (missingCol == g.cols) return; // no missing markers

  // 2. move non-missing sites to overwrite missing sites
  for (int col = missingCol + 1;
       col < g.cols; ++col) {
    if (hasMissingMarker(g, col)) {
      ++col;
    }

    for (int r = 0; r < g.rows; ++r ){
      g[r][missingCol] = g[r][col];
    }
    missingCol ++;
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
  int monoCol  = 0;
  for (monoCol = 0; monoCol < g.cols; ++monoCol) {
    if (isMonomorphicMarker(g, monoCol)) {
      break;
    }
  }
  if (monoCol == g.cols) return; // no monomorphic site

  // 2. move polymorphic sites to overwrite monomorphic sites
  for (int col = monoCol + 1; col < g.cols; ++col) {
    if (isMonomorphicMarker(g, col)) {
      continue;
    }

    // move g[, col] to g[, monoCol]
    for (int r = 0; r < g.rows; ++r) {
      g[r][monoCol] = g[r][col];
    }
    monoCol ++;
  }
  g.Dimension(g.rows, monoCol);
}

DataConsolidator::DataConsolidator()
    : strategy(DataConsolidator::UNINITIALIZED),
      phenotypeUpdated(true),
      covariateUpdated(true),
      kinshipForAuto(NULL),
      kinshipUForAuto(NULL),
      kinshipSForAuto(NULL),
      kinshipLoadedForAuto(false),
      kinshipForX(NULL),
      kinshipUForX(NULL),
      kinshipSForX(NULL),
      kinshipLoadedForX(false),
      sex(NULL),
      parRegion(NULL) {}

DataConsolidator::~DataConsolidator() {
  if (kinshipForAuto) delete kinshipForAuto;
  if (kinshipUForAuto) delete kinshipUForAuto;
  if (kinshipSForAuto) delete kinshipSForAuto;
  if (kinshipForX) delete kinshipForX;
  if (kinshipUForX) delete kinshipUForX;
  if (kinshipSForX) delete kinshipSForX;
}

/**
 * Load kinship file @param fn, load samples given by @param names
 * Store results to @param pKinship, @param pKinshipU, @param pKinshipS
 */
int DataConsolidator::loadKinshipFile(const std::string& fn,
                                      const std::vector<std::string>& names,
                                      EigenMatrix** pKinship,
                                      EigenMatrix** pKinshipU,
                                      EigenMatrix** pKinshipS,
                                      bool* pKinshipLoaded) {
  if (fn.empty()) {
    return -1;
  }
  if (!isUnique(names)) {
    return -1;
  }

  EigenMatrix*& kinship = *pKinship;
  EigenMatrix*& kinshipU = *pKinshipU;
  EigenMatrix*& kinshipS = *pKinshipS;
  bool& kinshipLoaded = *pKinshipLoaded;

  if (!kinship) {
    kinship = new EigenMatrix;
    kinshipU = new EigenMatrix;
    kinshipS = new EigenMatrix;
  }
  if (fn == "IDENTITY" || fn == "UNRELATED") {
    kinship->mat.resize(names.size(), names.size());
    kinship->mat.setZero();
    kinship->mat.diagonal().setOnes();
    kinship->mat.diagonal() *= 0.5;
    return 0;
  }

  LineReader lr(fn);
  int lineNo = 0;
  int fieldLen = 0;
  std::vector<std::string> fd;
  std::vector<int> columnToExtract;
  std::vector<std::string> header;  // kinship header line
  Eigen::MatrixXf& mat = kinship->mat;
  std::map<std::string, int> nameMap;
  makeMap(names, &nameMap);

  while (lr.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    if (lineNo == 1) {  // check header
      header = fd;
      fieldLen = fd.size();
      if (fieldLen < 2) {
        logger->error(
            "Insufficient column number (<2) in the first line of kinsihp "
            "file!");
        return -1;
      };
      if (tolower(fd[0]) != "fid" || tolower(fd[1]) != "iid") {
        logger->error("Kinship file header should begin with \"FID IID\"!");
        return -1;
      }
      std::map<std::string, int> headerMap;
      makeMap(fd, &headerMap);
      if (fd.size() != headerMap.size()) {
        logger->error("Kinship file have duplicated header!");
        return -1;
      }
      for (size_t i = 0; i < names.size(); ++i) {
        if (headerMap.count(names[i]) == 0) {
          logger->error(
              "The PID [ %s ] you specified cannot be found from kinship file!",
              names[i].c_str());
          return -1;
        }
        columnToExtract.push_back(headerMap[names[i]]);
      }
      mat.resize(columnToExtract.size(), columnToExtract.size());
      continue;
    }
    // body lines
    if ((int)fd.size() != fieldLen) {
      logger->error(
          "Inconsistent column number [ %zu ] (used to be [ %d ])in kinship "
          "file line [ %d ] - skip this file!",
          fd.size(), fieldLen, lineNo);
      return -1;
    }
    if (fd[1] != header[lineNo]) {  // PID in line i (1-based) should matched
                                    // the (i+1) (1-based) column in the header
      logger->error(
          "Inconsistent PID names in kinship file line [ %d ], file corrupted!",
          lineNo);
      return -1;
    }
    int row;
    if (nameMap.find(fd[1]) == nameMap.end()) {
      continue;
    }
    row = nameMap[fd[1]];
    for (size_t i = 0; i < columnToExtract.size(); ++i) {
      double d;
      if (str2double(fd[columnToExtract[i]], &d)) {
        mat(row, i) = d;
      } else {
        // unable to read, then set it to zero
        mat(row, i) = 0.0;
      }
    }
  }
#ifdef DEBUG
  std::string tmp = fn;
  tmp += ".tmp";
  std::ofstream ofs(tmp.c_str(), std::ofstream::out);
  ofs << mat;
  ofs.close();
#endif

  // fprintf(stderr, "Kinship matrix [ %d x %d ] loaded", (int)mat.rows(),
  // (int)mat.cols());
  kinshipLoaded = true;
  return 0;
}

/**
 * Decompose autosomal kinship matrix and release its memory
 * K = U * S * U'
 * @return 0 if success
 */
int DataConsolidator::decomposeKinshipForAuto() {
  if (!this->kinshipForAuto) return -1;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(this->kinshipForAuto->mat);
  if (es.info() == Eigen::Success) {
    (this->kinshipUForAuto->mat) = es.eigenvectors();
    (this->kinshipSForAuto->mat) = es.eigenvalues();

    if (this->kinshipForAuto) {
      delete this->kinshipForAuto;
      this->kinshipForAuto = NULL;
    }
    return 0;
  }
  return -1;
}

int DataConsolidator::decomposeKinshipForX() {
  if (!this->kinshipForX) return -1;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(this->kinshipForX->mat);
  if (es.info() == Eigen::Success) {
    (this->kinshipUForX->mat) = es.eigenvectors();
    (this->kinshipSForX->mat) = es.eigenvalues();

    if (this->kinshipForX) {
      delete this->kinshipForX;
      this->kinshipForX = NULL;
    }
    return 0;
  }
  return -1;
}

const EigenMatrix* DataConsolidator::getKinshipForAuto() const {
  if (!this->kinshipForAuto) return NULL;
  return this->kinshipForAuto;
}
const EigenMatrix* DataConsolidator::getKinshipUForAuto() const {
  if (!this->kinshipUForAuto) return NULL;
  return this->kinshipUForAuto;
}
const EigenMatrix* DataConsolidator::getKinshipSForAuto() const {
  if (!this->kinshipSForAuto) return NULL;
  return this->kinshipSForAuto;
}

const EigenMatrix* DataConsolidator::getKinshipForX() const {
  if (!this->kinshipForX) return NULL;
  return this->kinshipForX;
}
const EigenMatrix* DataConsolidator::getKinshipUForX() const {
  if (!this->kinshipUForX) return NULL;
  return this->kinshipUForX;
}
const EigenMatrix* DataConsolidator::getKinshipSForX() const {
  if (!this->kinshipSForX) return NULL;
  return this->kinshipSForX;
}

#if 0
  if (!this->kinshipLoadedForX) {
    this->kinshipForAutoAsKinshipForX = true;
    this->kinshipForAuto = this->kinshipForX;
    this->kinshipUForAuto = this->kinshipUForX;
    this->kinshipSForAuto = this->kinshipSForX;
    this->kinshipLoadedForX = true;
  }
#endif
