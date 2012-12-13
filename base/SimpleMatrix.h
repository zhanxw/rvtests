#ifndef _SIMPLEMATRIX_H_
#define _SIMPLEMATRIX_H_

#include <vector>
#include <string>
#include <stdio.h>

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
    for (unsigned int i = 0; i < d.size(); i++) {
      mat[i].push_back(d[i]);
    }
    return 0;
  }
  void clear() {
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
private:
  std::vector<std::string> rowName;
  std::vector<std::string> colName;
  std::vector< std::vector<double> > mat;
};

#endif /* _SIMPLEMATRIX_H_ */
