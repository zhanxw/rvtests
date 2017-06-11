#include "KinshipHolder.h"

#include <string>
#include <vector>

#include "third/eigen/Eigen/Core"
#include "third/eigen/Eigen/Eigenvalues"

#include "base/CommonFunction.h"
#include "base/IO.h"
#include "base/Logger.h"
#include "regression/EigenMatrix.h"
#define _USE_CXX11  // use C++11 timer
#include "base/SimpleTimer.h"
#include "base/TypeConversion.h"

#include "third/rapidjson/include/rapidjson/internal/dtoa.h"  // for

#undef DEBUG
//#define DEBUG
#ifdef DEBUG
#include <fstream>
#include <iostream>
#endif

extern Logger* logger;

Eigen::MatrixXf centerMatrix(const Eigen::MatrixXf& x) {
  return x.rowwise() - x.colwise().mean();
}

float corr(Eigen::MatrixXf& v1, Eigen::MatrixXf& v2) {
  Eigen::MatrixXf v1c = centerMatrix(v1);
  Eigen::MatrixXf v2c = centerMatrix(v2);

  float v1Norm = v1c.col(0).norm();
  float v2Norm = v2c.col(0).norm();
  if (v1Norm == 0.0 || v2Norm == 0.0) {
    return 0.0;
  }
  float v1v2 = v1c.col(0).dot(v2c.col(0));
  return v1v2 / (v1Norm * v2Norm);
}

KinshipHolder::KinshipHolder() {
  this->matK = new EigenMatrix;
  this->matU = new EigenMatrix;
  this->matS = new EigenMatrix;
  this->pSample = NULL;
  this->loaded = false;
}

KinshipHolder::~KinshipHolder() {
  release(&this->matK);
  release(&this->matU);
  release(&this->matS);
}

void KinshipHolder::release(EigenMatrix** m) {
  if (m && *m) {
    delete (*m);
    *m = NULL;
  }
}

int KinshipHolder::setSample(const std::vector<std::string>& names) {
  this->pSample = &names;
  if (!isUnique(names)) {
    return -1;
  }
  return 0;
}

int KinshipHolder::setFile(const std::string& fileName) {
  this->fileName = fileName;

  // kinship file may not exists
  // e.g. allow "IDENTITY" or "UNRELATED" as the file name and represents an
  // identity matrix

  return 0;
}

int KinshipHolder::setEigenFile(const std::string& fileName) {
  this->eigenFileName = fileName;

  // eigen file may not exists
  // if not, we will write it with eigen-decomposition results.

  return 0;
}

int KinshipHolder::load() {
  AccurateTimer timer;
  // deal with special matrix
  if (isSpecialFileName()) {
    loadDecomposedIdentityKinship();
    logger->info(
        "DONE: Loaded special kinship [ %s ] successfully in [ %.1f ] seconds.",
        this->fileName.c_str(), timer.stop());
    this->loaded = true;
    return 0;
  }

  int ret = 0;
  // load kinship
  ret = loadK();
  if (ret) {
    logger->error("Failed to load kinship file [ %s ]", this->fileName.c_str());
    return -1;
  }
  logger->info(
      "DONE: Loaded kinship file [ %s ] successfully in [ %.1f ] seconds.",
      this->fileName.c_str(), timer.stop());

  // load kinship eigen file
  bool isDecomposed = false;
  if (this->eigenFileName.size() && fileExists(this->eigenFileName)) {
    timer.start();
    ret = loadDecomposed();
    if (!ret) {
      logger->info(
          "DONE: Loaded spectral decomposition result file [ %s ] succeeded "
          "in [ %.1f ] seconds.",
          this->eigenFileName.c_str(), timer.stop());
      isDecomposed = true;
    } else {
      logger->warn(
          "Failed to load spectral decomposition results of the kinship "
          "matrix!");
    }
  }
  if (!isDecomposed) {
    timer.start();
    ret = decompose();

    if (ret) {
      logger->error("Failed to decompose kinship matrix");
      return -1;
    }
    logger->info(
        "DONE: Spectral decomposition of the kinship matrix succeeded in [ "
        "%.1f ] seconds.",
        timer.stop());
    this->loaded = true;
  }

  if (!fileExists(this->eigenFileName) && this->eigenFileName.size()) {
    timer.start();
    ret = saveDecomposed();
    if (ret) {
      logger->error("Cannot store spectral decomposition results [ %s ]",
                    this->eigenFileName.c_str());
    }
  }
  return 0;
}

/**
 * Load kinship file and store the kinship matrix internally
 */
int KinshipHolder::loadK() {
  // deal with special matrix
  if (isSpecialFileName()) {
    loadIdentityKinship();
    return 0;
  }

  // check kinship file existance
  if (!fileExists(this->fileName)) {
    logger->warn("Cannot open kinship file [ %s ]", this->fileName.c_str());
    return -1;
  }

  LineReader lr(this->fileName);
  int lineNo = 0;
  int fieldLen = 0;
  std::vector<std::string> fd;
  std::vector<int> columnToExtract;
  std::vector<std::string> header;  // kinship header line
  Eigen::MatrixXf& mat = this->matK->mat;
  const std::vector<std::string>& names = *this->pSample;
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
          logger->error("Cannot find sample [ %s ] from the kinship file!",
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
          "Inconsistent column number [ %zu ] (used to be [ %d ]) in kinship "
          "file line [ %d ] - skip this file!",
          fd.size(), fieldLen, lineNo);
      // todo: this error may be false alarm - as the last line may be empty
      continue;
      // exit(1);
    }
    if (fd[1] != header[lineNo]) {  // PID in line i (1-based) should matched
                                    // the (i+1) (1-based) column in the header
      logger->error(
          "Inconsistent IID names in kinship file line [ %d ], file corrupted!",
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
  this->loaded = true;
  return 0;
}

int KinshipHolder::decompose() {
  // eigen decomposition
  if (!this->matK) {
    fprintf(stderr,
            "%s:%d cannot dereference and decompose a null-pointed matrix!\n",
            __FILE__, __LINE__);
    return -1;
  }
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(this->matK->mat);
  if (es.info() == Eigen::Success) {
    (this->matU->mat) = es.eigenvectors();
    (this->matS->mat) = es.eigenvalues();

    if (this->matK) {
      delete this->matK;
      this->matK = NULL;
    }
    return 0;
  }
  return -1;
}

int KinshipHolder::loadDecomposed() {
  LineReader lr(this->eigenFileName);
  int lineNo = 0;
  int fieldLen = 0;
  std::vector<std::string> fd;
  std::vector<int> columnToExtract;
  std::vector<std::string> header;  // header line of the kinship eigen file
  Eigen::MatrixXf& matK = this->matK->mat;
  Eigen::MatrixXf& matS = this->matS->mat;
  Eigen::MatrixXf& matU = this->matU->mat;
  const std::vector<std::string>& names = *this->pSample;
  const int NumSample = (int)names.size();
  std::map<std::string, int> nameMap;
  makeMap(names, &nameMap);
  std::map<std::string, int> headerMap;

  while (lr.readLineBySep(&fd, "\t ")) {
    ++lineNo;
    if (lineNo == 1) {  // check header
      header = fd;
      fieldLen = fd.size();
      if (fieldLen < 3) {  // at least three columns: IID, Lambda, U1
        logger->error(
            "Insufficient column number (<3) in the first line of kinsihp "
            "file!");
        return -1;
      };
      for (size_t i = 0; i != fd.size(); ++i) {
        fd[i] = tolower(fd[i]);
      }
      makeMap(fd, &headerMap);
      if (fd.size() != headerMap.size()) {
        logger->error("Kinship file have duplicated headers!");
        return -1;
      }

      // check IID, Lambda, U1, U2, ... U(N) where (N) is the sample size
      if (headerMap.count("iid") == 0) {
        logger->error("Missing 'IID' column!");
        return -1;
      }
      columnToExtract.push_back(headerMap["iid"]);

      if (headerMap.count("lambda") == 0) {
        logger->error("Missing 'Lambda' column!");
        return -1;
      }
      columnToExtract.push_back(headerMap["lambda"]);

      std::string s;
      for (int i = 0; i < NumSample; ++i) {
        s = "u";
        s += toString(i + 1);
        if (headerMap.count(s) == 0) {
          logger->warn(
              "Missing [ %s ] column in the header line of file [ %s ] when we "
              "are analyzing [ %d ] "
              "samples!",
              s.c_str(), this->eigenFileName.c_str(), NumSample);
          return -1;
        }
        columnToExtract.push_back(headerMap[s]);
      }
      s = "u";
      s += toString(NumSample + 1);
      if (headerMap.count(s) != 0) {
        logger->error("Unexpected column '%s'!", s.c_str());
        return -1;
      }

      matS.resize(NumSample, 1);
      matU.resize(NumSample, NumSample);
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

    const int iidColumn = columnToExtract[0];
    const std::string& iid = fd[iidColumn];
    if (nameMap.count(iid) == 0) {
      logger->error("Unexpected sample [ %s ]!", iid.c_str());
      return -1;
    }
    const int row = nameMap[iid];

    const int lambdaColumn = columnToExtract[1];
    double temp = 0.0;
    if (!str2double(fd[lambdaColumn], &temp)) {
      logger->warn("Invalid numeric value [ %s ] treated as zero!",
                   fd[lambdaColumn].c_str());
    }
    matS(lineNo - 2, 0) = temp;

    for (int i = 0; i < NumSample; ++i) {
      int uColumn = columnToExtract[i + 2];
      if (!str2double(fd[uColumn], &temp)) {
        logger->warn("Invalid numeric value [ %s ] treated as zero!",
                     fd[lambdaColumn].c_str());
      }
      matU(row, i) = temp;
    }
  }

  // verify eigen decomposition results make senses
  // check largest eigen vector and eigen value
  Eigen::MatrixXf v1 = matK * matU.col(0);
  Eigen::MatrixXf v2 = matS(0, 0) * matU.col(0);
  if (matS(0, 0) > 0.5 && v1.col(0).norm() > .5 && v2.col(0).norm() > 0.5 &&
      corr(v1, v2) < 0.8) {
    logger->warn("Cannot verify spectral decompose results!");
    return -1;
  }

  // check the min(10, NumSample) random eigen vector and eigen value
  int randomCol = 10;
  if (randomCol > NumSample - 1) {
    randomCol = NumSample - 1;
  }
  v1 = matK * matU.col(randomCol);
  v2 = matS(randomCol, 0) * matU.col(randomCol);
  if (matS(randomCol, 0) > 0.5 && v1.col(0).norm() > 0.5 &&
      v2.col(0).norm() > 0.5 && corr(v1, v2) < 0.8) {
    logger->warn("Cannot verify spectral decompose results!");
    return -1;
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

  if (this->matK) {
    delete this->matK;
    this->matK = NULL;
  }

  return 0;
}

int KinshipHolder::saveDecomposed() {
  char buffer[1024];
  FileWriter fw(this->eigenFileName.c_str());
  const std::vector<std::string>& names = *this->pSample;
  const int nSample = (int)names.size();

  fw.write("IID\tLambda");
  for (int i = 0; i < nSample; ++i) {
    fw.printf("\tU%d", i + 1);
  }
  fw.write("\n");

  for (int i = 0; i < nSample; ++i) {
    fw.write(names[i].c_str());
    fw.write("\t");
    char* p = rapidjson::internal::dtoa(matS->mat(nSample - 1 - i), buffer);
    *p = '\0';
    fw.write(buffer);
    for (int j = 0; j < nSample; ++j) {
      char* p =
          rapidjson::internal::dtoa(matU->mat(i, nSample - 1 - j), buffer);
      *p = '\0';
      fw.write("\t");
      fw.write(buffer);
    }
    fw.write("\n");
  }
  return 0;
}

bool KinshipHolder::isSpecialFileName() {
  if (fileName == "IDENTITY" || fileName == "UNRELATED") {
    return true;
  }
  return false;
}

int KinshipHolder::loadIdentityKinship() {
  const std::vector<std::string>& names = *this->pSample;
  assert(isSpecialFileName());

  this->matK->mat.resize(names.size(), names.size());
  this->matK->mat.setZero();
  this->matK->mat.diagonal().setOnes();

  return 0;
}

int KinshipHolder::loadDecomposedIdentityKinship() {
  const std::vector<std::string>& names = *this->pSample;
  assert(isSpecialFileName());

  // this->matK->mat.resize(names.size(), names.size());
  // this->matK->mat.setZero();
  // this->matK->mat.diagonal().setOnes();
  // this->matK->mat.diagonal() *= 0.5;

  const int N = names.size();
  this->matU->mat = Eigen::MatrixXf::Identity(N, N);
  this->matS->mat = Eigen::VectorXf::Ones(N);

  return 0;
}
