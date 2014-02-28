#ifndef _SIMPLEMATRIX_H_
#define _SIMPLEMATRIX_H_

#include <vector>
#include <string>
#include <stdio.h>

/**
 * This matrix class is for convenient store matrix class.
 */
class SimpleMatrix{
 public:
  SimpleMatrix() {};
  SimpleMatrix(int nr, int nc) {
    this->resize(nr, nc);
  }
  /* const static int COLUMN_HEADER = 0x1; */
  /* const static int ROW_HEADER = 0x2; */
  int readFile(const char* f);
  int writeFile(const char* f);
  std::vector<double>& operator[] (int i) {
    return mat[i];
  }
  const std::vector<double>& operator[] (int i) const {
    return mat[i];
  }
  void resize(int nrow, int ncol) {
    if (nrow < 0 || ncol < 0) return;
    mat.resize(nrow);
    for (int i = 0; i < nrow; i++)
      mat[i].resize(ncol);
  };
  int appendRow(const std::vector<double>& d) {
    if (mat.size() && mat[0].size() != d.size()) {
      fprintf(stderr, "Row size does not match!\n");
      return -1;
    }
    mat.push_back(d);
    return 0;
  };
  int appendCol(const std::vector<double>& d) {
    if (mat.size()) {
      // no empty
      if (mat.size() != d.size()) {
        fprintf(stderr, "Col size does not match!\n");
        return -1;
      }
    }

    mat.resize(d.size());
    for (size_t i = 0; i < d.size(); i++) {
      mat[i].push_back(d[i]);
    }
    return 0;
  }
  int deleteRow(int i) {
    if (i < 0 || (size_t)i > mat.size()) return -1;
    mat.erase(mat.begin() + i);
    if (!rowName.empty())
      rowName.erase(rowName.begin() + i);
    return 0;
  }
  int deleteCol(int i) {
    if (i < 0 || i > ncol()) return -1;
    for (size_t nr = 0; nr < mat.size(); ++nr) {
      mat[nr].erase(mat[nr].begin() + i);
    }
    if (!colName.empty()){
      colName.erase(colName.begin() + i);
    }
    return 0;
  }
  void clear() {
    mat.clear();
  }
  void zero() {
    for (unsigned int i = 0; i < mat.size(); i++){
      for (unsigned int j = 0; j < mat[i].size(); j++){
        mat[i][j] = 0.0;
      }
    }
  }
  int nrow() const {
    return mat.size();
  };
  int ncol() const {
    if (mat.size() == 0)
      return 0;
    return mat[0].size();
  };
  const std::vector<std::string>& getRowName() const{ return this->rowName;};
  const std::vector<std::string>& getColName() const{ return this->colName;};
  int setRowName(const int idx, const std::string& s){
    if (idx < 0 || idx >= nrow()) return -1;
    rowName.resize(nrow());
    rowName[idx] = s;
    return 0;
  }
  int setColName(const int idx, const std::string& s){
    if (idx < 0 || idx >= ncol()) return -1;
    colName.resize(ncol());
    colName[idx] = s;
    return 0;
  }
 private:
  std::vector<std::string> rowName;
  std::vector<std::string> colName;
  std::vector< std::vector<double> > mat;
};

#endif /* _SIMPLEMATRIX_H_ */
