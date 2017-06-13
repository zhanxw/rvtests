#include "SimpleMatrix.h"

#include <algorithm>
#include <cmath>

#include "CommonFunction.h"
#include "IO.h"
#include "Indexer.h"
#include "TypeConversion.h"

void SimpleMatrix::resize(int nr, int nc) {
  if (nr < 0 || nc < 0) return;

  const int onr = nrow();
  const int onc = ncol();

  mat.resize(nr);
  for (int i = 0; i < nr; i++) mat[i].resize(nc);
  rowName.resize(nr);
  colName.resize(nc);

  for (int i = onr; i < nr; ++i) {
    rowName[i] = ("R");
    rowName[i] += toString(i + 1);
  }
  for (int i = onc; i < nc; ++i) {
    colName[i] = "C";
    colName[i] += toString(i + 1);
  }
}

/**
 * @return 0: success
 *
 * will read row names and column names.
 */
int SimpleMatrix::readFile(const char* f) {
  this->clear();
  LineReader lr(f);
  std::vector<std::string> fd;
  std::vector<double> d;
  int lineNo = 0;
  while (lr.readLineBySep(&fd, " \t")) {
    lineNo++;
    if (lineNo == 1) {
      colName.resize(fd.size());
      std::copy(fd.begin() + 1, fd.end(), colName.begin());
    } else {
      d.resize(fd.size() - 1);
      for (size_t i = 1; i != fd.size(); i++) d[i - 1] = atof(fd[i]);
      if (this->nrow() > 0 && this->ncol() + 1 != (int)fd.size()) {
        fprintf(stderr,
                "At line %d, column width is not consistent (expected: %d, "
                "actual: %d )!\n",
                lineNo, this->ncol(), (int)fd.size());
        continue;
      }
      rowName.push_back(fd[0]);
      this->mat.push_back(d);
    }
  }
  return 0;
}

/**
 * @return 0: success
 *
 * will write row names and column names.
 */
int SimpleMatrix::writeFile(const char* f) const {
  FileWriter fw(f);
  const int nr = this->nrow();
  const int nc = this->ncol();
  fw.printf("R%dC%d", nr, nc);
  for (int j = 0; j < nc; ++j) {
    fw.write("\t");
    fw.write(colName[j]);
  }
  fw.write("\n");

  for (int i = 0; i < nr; i++) {
    fw.write(rowName[i]);
    for (int j = 0; j < nc; j++) {
      fw.write("\t");
      fw.printf("%f", mat[i][j]);
    }
    fw.write("\n");
  }
  return 0;
}

template <typename T>
int SimpleMatrix::appendRow(const std::vector<T>& d, const std::string& label) {
  if (mat.size() && mat[0].size() != d.size()) {
    fprintf(stderr, "Column size does not match!\n");
    return -1;
  }
  mat.push_back(d);

  if (label.empty()) {
    rowName.push_back("R");
    rowName[rowName.size() - 1] += toString(rowName.size());
  } else {
    rowName.push_back(label);
  }

  if (!colName.size()) {
    resetColName();
  }
  return 0;
}
template int SimpleMatrix::appendRow(const std::vector<double>& d,
                                     const std::string& label);

template <typename T>
int SimpleMatrix::appendCol(const std::vector<T>& d, const std::string& label) {
  if (mat.size()) {
    // no empty
    if (mat.size() != d.size()) {
      fprintf(stderr, "Row size does not match!\n");
      return -1;
    }
  } else {
    // empty
    mat.resize(d.size());
    resetRowName();
  }

  for (size_t i = 0; i < d.size(); i++) {
    mat[i].push_back(d[i]);
  }

  if (label.empty()) {
    colName.push_back("C");
    colName[colName.size() - 1] += toString(colName.size());
  } else {
    colName.push_back(label);
  }

  return 0;
}
template int SimpleMatrix::appendCol(const std::vector<double>& d,
                                     const std::string& label);
template int SimpleMatrix::appendCol(const std::vector<int>& d,
                                     const std::string& label);

int SimpleMatrix::deleteRow(int i) {
  if (i < 0 || (size_t)i > mat.size()) return -1;
  mat.erase(mat.begin() + i);
  if (!rowName.empty()) rowName.erase(rowName.begin() + i);
  return 0;
}
int SimpleMatrix::deleteCol(int i) {
  if (i < 0 || i > ncol()) return -1;
  for (size_t nr = 0; nr < mat.size(); ++nr) {
    mat[nr].erase(mat[nr].begin() + i);
  }
  if (!colName.empty()) {
    colName.erase(colName.begin() + i);
  }
  return 0;
}

int SimpleMatrix::setRowName(const int idx, const std::string& s) {
  if (idx < 0 || idx >= nrow()) return -1;
  // rowName.resize(nrow());
  rowName[idx] = s;
  return 0;
}
int SimpleMatrix::setColName(const int idx, const std::string& s) {
  if (idx < 0 || idx >= ncol()) return -1;
  // colName.resize(ncol());
  colName[idx] = s;
  return 0;
}

int SimpleMatrix::setRowName(const std::vector<std::string>& name) {
  if (rowName.size() == name.size()) {
    rowName = name;
    return 0;
  }
  return -1;
}
int SimpleMatrix::setColName(const std::vector<std::string>& name) {
  if (colName.size() == name.size()) {
    colName = name;
    return 0;
  }
  return -1;
}

int SimpleMatrix::dropRow(const std::set<std::string>& rowNameSet) {
  if (rowNameSet.empty()) {
    return 0;
  }

  const int nr = nrow();
  const int nc = ncol();
  int idx = 0;
  for (int i = 0; i < nr; ++i) {
    if (rowNameSet.count(rowName[i])) {  // need to drop this row
      continue;
    }
    assignRow(idx, i);
    idx++;
  }
  resize(idx, nc);
  return 0;
}

int SimpleMatrix::dropRow(const std::vector<std::string>& name) {
  std::set<std::string> s(name.begin(), name.end());
  return dropRow(s);
}

int SimpleMatrix::dropRow(const std::vector<int>& index) {
  removeByIndex(index, &mat);
  removeByIndex(index, &rowName);
  return 0;
}

int SimpleMatrix::addRow(const std::vector<std::string>& newRowName,
                         double value) {
  const int nc = ncol();
  std::vector<double> fd(nc, value);
  extend(newRowName, &rowName);
  for (size_t i = 0; i != newRowName.size(); ++i) {
    mat.push_back(fd);
  }
  return 0;
}

int SimpleMatrix::keepRow(const std::vector<std::string>& name) {
  std::set<std::string> s(name.begin(), name.end());
  std::vector<int> indexToRemove;
  const int nr = nrow();
  for (int i = 0; i < nr; ++i) {
    if (!s.count(rowName[i])) {
      indexToRemove.push_back(i);
    }
  }
  removeByIndex(indexToRemove, &rowName);
  removeByIndex(indexToRemove, &mat);
  return 0;
}

int SimpleMatrix::keepCol(const std::vector<std::string>& name) {
  const int nc = ncol();
  const int nr = nrow();

  std::vector<int> indexToRemove;
  std::set<std::string> keep(name.begin(), name.end());

  for (int i = 0; i < nc; ++i) {
    if (!keep.count(colName[i])) {
      indexToRemove.push_back(i);
    }
  }

  removeByIndex(indexToRemove, &this->colName);
  for (int i = 0; i < nr; ++i) {
    removeByIndex(indexToRemove, &this->mat[i]);
  }
  return 0;
}

int SimpleMatrix::assignRow(const int to, const int from) {
  if (to < 0 || from < 0 || to >= (int)rowName.size() ||
      from >= (int)rowName.size()) {
    return -1;
  }
  mat[to] = mat[from];
  rowName[to] = rowName[from];
  return 0;
}

void SimpleMatrix::extractCol(int col, std::vector<double>* v) const {
  if (!v) return;
  if (col >= ncol()) return;
  int nr = nrow();
  v->resize(nr);
  for (int i = 0; i < nr; ++i) {
    (*v)[i] = mat[i][col];
  }
}

std::vector<double> SimpleMatrix::extractCol(int col) const {
  std::vector<double> v;
  extractCol(col, &v);
  return v;
}

int SimpleMatrix::setCol(int col, const std::vector<double>& v) {
  int nr = nrow();
  assert(nr == (int)v.size());
  assert(col < (int)ncol());

  for (int i = 0; i < nr; ++i) {
    mat[i][col] = v[i];
  }

  return 0;
}

void SimpleMatrix::resetRowName() {
  const int nrow = this->nrow();
  rowName.resize(nrow);
  for (int i = 0; i < nrow; ++i) {
    rowName[i] = ("R");
    rowName[i] += toString(i + 1);
  }
}
void SimpleMatrix::resetColName() {
  const int ncol = this->ncol();
  colName.resize(ncol);
  for (int i = 0; i < ncol; ++i) {
    colName[i] = "C";
    colName[i] += toString(i + 1);
  }
}

std::vector<int> SimpleMatrix::allMissingRows() const {
  std::vector<int> ret;
  const int nr = nrow();
  const int nc = ncol();
  bool allMissing = true;
  for (int i = 0; i < nr; ++i) {
    allMissing = true;
    for (int j = 0; j < nc; ++j) {
      if (!std::isnan(mat[i][j])) {
        allMissing = false;
        break;
      }
    }
    if (allMissing) {
      ret.push_back(i);
    }
  }
  return ret;
}

int SimpleMatrix::reorderRow(const std::vector<int>& indice) {
  std::vector<std::vector<double> > newMat;
  std::vector<std::string> newRowName;
  const size_t n = indice.size();
  const int nr = nrow();
  for (size_t i = 0; i != n; ++i) {
    if (indice[i] < 0 || indice[i] >= nr) {
      continue;
    }
    newMat.push_back(mat[indice[i]]);
    newRowName.push_back(rowName[indice[i]]);
  }

  std::swap(mat, newMat);
  std::swap(rowName, newRowName);

  return 0;
}

int SimpleMatrix::reorderRow(const std::vector<std::string>& indice) {
  std::vector<int> ind;
  Indexer indexer(rowName);
  indexer.translate(indice, &ind);
  return reorderRow(ind);
}

int SimpleMatrix::imputeMissingToMeanByCol() {
  const int nr = nrow();
  const int nc = ncol();
  std::vector<double> d(nc, 0.0);
  std::vector<int> n(nc, 0);
  std::vector<int> hasMissing(nc, false);

  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      if (std::isfinite(mat[i][j])) {
        d[j] += mat[i][j];
        n[j]++;
      } else {
        hasMissing[j] = true;
      }
    }
  }

  double avg;
  for (int j = 0; j < nc; ++j) {
    if (!hasMissing[j]) continue;
    avg = 0.0;
    if (n[j] != 0) {
      avg = d[j] / n[j];
    }
    for (int i = 0; i < nr; ++i) {
      if (!std::isfinite(mat[i][j])) {
        mat[i][j] = avg;
      }
    }
  }

  return 0;
}
