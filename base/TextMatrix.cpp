#include "TextMatrix.h"

#include "CommonFunction.h"
#include "IO.h"
#include "SimpleMatrix.h"
#include "TypeConversion.h"

int TextMatrix::readFile(const std::string& fn, int flag) {
  this->clear();
  LineReader lr(fn);
  std::vector<std::string> fd;
  std::vector<std::string> d;
  int lineNo = 0;
  while (lr.readLineBySep(&fd, " \t")) {
    lineNo++;
    if (lineNo == 1 && flag & TextMatrix::HAS_HEADER) {
      colName.resize(fd.size());
      if (flag & TextMatrix::HAS_ROWNAME) {
        std::copy(fd.begin() + 1, fd.end(), colName.begin());
      } else {
        colName = fd;
      }
    } else {
      if (flag & TextMatrix::HAS_ROWNAME) {
        d.resize(fd.size() - 1);
        for (size_t i = 1; i != fd.size(); i++) d[i - 1] = fd[i];
        if (this->nrow() > 0 && this->ncol() + 1 != (int)fd.size()) {
          fprintf(stderr,
                  "At line %d, column width is not consistent (expected: %d, "
                  "actual: %d )!\n",
                  lineNo, this->ncol(), (int)fd.size());
          continue;
        }
        rowName.push_back(fd[0]);
        this->mat.push_back(d);
      } else {
        if (lineNo != 1 && fd.size() != this->colName.size()) {
          fprintf(
              stderr,
              "At line %d, the number of columns does not equal to [ %d ]!\n",
              lineNo, (int)colName.size());
          return -1;
        }
        this->mat.push_back(fd);
        this->rowName.push_back("R");
        this->rowName.back() += toString(rowName.size());
      }
    }
  }
  assert(this->rowName.size() == (size_t) this->nrow());
  assert(this->colName.size() == (size_t) this->ncol());
  return 0;
}

int TextMatrix::writeFile(const std::string& fn, int flag) const {
  FileWriter fw(fn);
  const int NR = nrow();
  const int NC = ncol();
  if (flag & TextMatrix::OUTPUT_HEADER) {
    if (flag & TextMatrix::OUTPUT_ROWNAME) {
      fw.printf("R%dC%d", NR, NC);
    }
    for (size_t i = 0; i != colName.size(); ++i) {
      fw.write("\t");
      fw.write(colName[i].c_str());
    }
    fw.write("\n");
  }
  for (int i = 0; i < NR; ++i) {
    if (flag & TextMatrix::OUTPUT_ROWNAME) {
      fw.write(rowName[i].c_str());
      for (int j = 0; j < NC; ++j) {
        fw.write("\t");
        fw.write(mat[i][j].c_str());
      }
    } else {  // no rowname
      for (int j = 0; j < NC; ++j) {
        if (j) {
          fw.write("\t");
        }
        fw.write(mat[i][j].c_str());
      }
    }
    fw.write("\n");
  }
  return 0;
}

const std::vector<std::string> TextMatrix::header() const { return colName; }

int TextMatrix::keepCol(const std::vector<std::string>& name) {
  const int nc = ncol();
  const int nr = nrow();

  std::vector<int> indexToRemove;
  std::set<std::string> keep;
  makeSet(name, &keep);

  for (int i = 0; i < nc; ++i) {
    if (keep.count(colName[i])) continue;
    indexToRemove.push_back(i);
  }

  removeByIndex(indexToRemove, &colName);
  for (int i = 0; i < nr; ++i) {
    removeByIndex(indexToRemove, &mat[i]);
  }
  assert(name.size() == (size_t)ncol());
  return 0;
}

int TextMatrix::setRowNameByCol(const std::string& name) {
  const int nc = ncol();
  const int nr = nrow();
  int idx = -1;
  for (int i = 0; i < nc; ++i) {
    if (colName[i] == name) {
      idx = i;
      break;
    }
  }
  if (idx < 0) return -1;
  for (int i = 0; i < nr; ++i) {
    rowName[i] = mat[i][idx];
  }
  return 0;
}

int TextMatrix::convert(SimpleMatrix* out) const {
  const int nc = ncol();
  const int nr = nrow();

  double d = 0;
  SimpleMatrix& m = *out;
  m.resize(nr, nc);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      if (str2double(mat[i][j], &d)) {
        m[i][j] = d;
      } else {
        m[i][j] = NAN;
      }
    }
  }

  m.setRowName(rowName);
  m.setColName(colName);

  return 0;
}

void TextMatrix::extractCol(int col, std::vector<std::string>* v) const {
  if (!v) return;
  if (col >= ncol()) return;
  int nr = nrow();
  v->resize(nr);
  for (int i = 0; i < nr; ++i) {
    (*v)[i] = mat[i][col];
  }
}

std::vector<std::string> TextMatrix::extractCol(int col) const {
  std::vector<std::string> v;
  extractCol(col, &v);
  return v;
}

void TextMatrix::extractCol(const std::string& col,
                            std::vector<std::string>* v) const {
  if (!v) return;
  int idx = -1;
  for (size_t i = 0; i != colName.size(); ++i) {
    if (colName[i] == col) {
      idx = i;
      break;
    }
  }
  if (idx < 0) return;

  int nr = nrow();
  v->resize(nr);
  for (int i = 0; i < nr; ++i) {
    (*v)[i] = mat[i][idx];
  }
}

std::vector<std::string> TextMatrix::extractCol(const std::string& col) const {
  std::vector<std::string> v;
  extractCol(col, &v);
  return v;
}
