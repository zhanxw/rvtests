#include "DataConsolidator.h"

#include "third/eigen/Eigen/Core"
#include "third/eigen/Eigen/Eigenvalues"

#include "CommonFunction.h"
#include "regression/EigenMatrix.h"
#include "IO.h"

DataConsolidator::DataConsolidator()
  :
  strategy(DataConsolidator::UNINITIALIZED),
  kinship(NULL), kinshipU(NULL), kinshipS(NULL), kinshipLoaded(false)
{
};
DataConsolidator::~DataConsolidator(){
  if (!kinship) delete kinship;
  if (!kinshipU) delete kinshipU;
  if (!kinshipS) delete kinshipS;
}
int DataConsolidator::loadKinshipFile(const std::string& fn, const std::vector<std::string>& names){
  if (fn.empty()) {
    return -1;
  }
  if (!isUnique(names)) {
    return -1;
  }
  if (!kinship) {
    kinship = new EigenMatrix;
    kinshipU = new EigenMatrix;
    kinshipS = new EigenMatrix;
  }
  if (fn == "IDENTITY" || fn == "UNRELATED") {
    kinship->mat.resize(names.size(), names.size());
    kinship->mat.setZero();
    kinship->mat.diagonal().setOnes();
    kinship->mat *= 0.5;
    return 0;
  }
  
  LineReader lr(fn);
  int lineNo = 0;
  int fieldLen = 0;
  std::vector<std::string> fd;
  std::vector<int> columnToExtract;
  std::vector<std::string> header; // kinship header line
  Eigen::MatrixXf& mat = this->kinship->mat;
  std::map<std::string, int> nameMap;
  makeMap(names, &nameMap);
  
  while (lr.readLineBySep(&fd, "\t ")){
    ++ lineNo;
    if (lineNo == 1) { //check header
      header = fd;
      fieldLen = fd.size();
      if (fieldLen < 2) {
        logger->error("Insufficient column number (<2) in the first line of kinsihp file!");
        return -1;
      };
      if (tolower(fd[0]) != "fid" ||
          tolower(fd[1]) != "iid") {
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
          logger->error("The PID [ %s ] you specified cannot be found from kinship file!", names[i].c_str());
          return -1;
        }
        columnToExtract.push_back(headerMap[names[i]]);
      }
      mat.resize(columnToExtract.size(), columnToExtract.size());
      continue;
    }
    // body lines
    if ((int)fd.size() != fieldLen) {
      logger->error("Inconsistent column number [ %zu ] (used to be [ %d ])in kinship file line [ %d ] - skip this file!", fd.size(), fieldLen, lineNo);
      return -1;
    }
    if (fd[1] != header[lineNo]) {// PID in line i (1-based) should matched the (i+1) (1-based) column in the header
      logger->error("Inconsistent PID names in kinship file line [ %d ], file corrupted!", lineNo);
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
  // fprintf(stderr, "Kinship matrix [ %d x %d ] loaded", (int)mat.rows(), (int)mat.cols());
  this->kinshipLoaded = true;
  return 0;
}
int DataConsolidator::decomposeKinship(){
  if (!this->kinship) return -1;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es(this->kinship->mat);
  if (es.info() == Eigen::Success) {
    (this->kinshipU->mat) = es.eigenvectors();
    (this->kinshipS->mat) = es.eigenvalues();

    if (this->kinship) {
      delete this->kinship;
      this->kinship = NULL;
    }
    return 0;  
  }
  return -1;
}

const EigenMatrix* DataConsolidator::getKinship() const{
  if (!this->kinship) return NULL;  
  return this->kinship;
}
const EigenMatrix* DataConsolidator::getKinshipU() const{
  if (!this->kinshipU) return NULL;  
  return this->kinshipU;
}
const EigenMatrix* DataConsolidator::getKinshipS() const{
    if (!this->kinshipS) return NULL;
  return this->kinshipS;
}
