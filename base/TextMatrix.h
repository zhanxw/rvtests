#ifndef TEXTMATRIX_H
#define TEXTMATRIX_H

#include <string>
#include <vector>

class SimpleMatrix;

class TextMatrix {
 public:
  /**
   * read in a text file
   */
  int readFile(const std::string& fn, int flag = 0);
  std::vector<std::string>& operator[](int i) { return mat[i]; }
  const std::vector<std::string> header() const;
  void clear() {
    mat.clear();
    rowName.clear();
    colName.clear();
  }
  int nrow() const { return mat.size(); };
  int ncol() const {
    if (mat.size() == 0) return 0;
    return mat[0].size();
  }
  const std::vector<std::string>& getRowName() const { return this->rowName; };
  const std::vector<std::string>& getColName() const { return this->colName; };

  /**
   * Only keep the columns whose names are in @param colName
   */
  int keepCol(const std::vector<std::string>& colName);

  int setRowNameByCol(const std::string& name);

  /**
   * Convert to SimpleMatrix @param m
   * For non-numeric values, a NAN value will be used.
   */
  int convert(SimpleMatrix* m) const;

  void extractCol(int col, std::vector<std::string>* v) const;
  std::vector<std::string> extractCol(int col) const;
  void extractCol(const std::string& col, std::vector<std::string>* v) const;
  std::vector<std::string> extractCol(const std::string& col) const;

 public:
  const static int HAS_HEADER = 1;
  const static int HAS_ROWNAME = 2;

 private:
  std::vector<std::string> rowName;
  std::vector<std::string> colName;
  std::vector<std::vector<std::string> > mat;
};

#endif /* TEXTMATRIX_H */
