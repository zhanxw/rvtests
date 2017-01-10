#ifndef _MATRIXIO_H_
#define _MATRIXIO_H_

#include "base/TypeConversion.h"

void LoadVector(const char* fn, Vector& v) {
  LineReader lr(fn);
  std::vector<std::string> s;
  int lineNo = 0;
  while (lr.readLineBySep(&s, " \t")) {
    lineNo++;
    v.Dimension(lineNo);
    // v[lineNo - 1] = atof(s[0].c_str());
    double d;
    if (str2double(s[0], &d)) {
      v[lineNo - 1] = d;
    } else {
      v[lineNo - 1] = NAN;
    }
  }
}

void LoadMatrix(const char* fn, Matrix& m) {
  LineReader lr(fn);
  std::vector<std::string> s;
  int lineNo = 0;
  while (lr.readLineBySep(&s, " \t")) {
    lineNo++;
    m.Dimension(lineNo, s.size());
    for (int j = 0; j < s.size(); j++) {
      // m[lineNo - 1][j] = atof(s[j].c_str());
      double d;
      if (str2double(s[j], &d)) {
        m[lineNo - 1][j] = d;
      } else {
        m[lineNo - 1][j] = NAN;
      }
    }
  }
}

void Print(const Vector& v) {
  for (int i = 0; i < v.Length(); i++) {
    if (i) {
      fprintf(stdout, "\t");
    }
    fprintf(stdout, "%.3f", v[i]);
  }
}

void Print(const Matrix& m) {
  for (int i = 0; i < m.rows; i++) {
    for (int j = 0; j < m.cols; j++) {
      if (j) {
        fprintf(stdout, "\t");
      }
      fprintf(stdout, "%.3f", m[i][j]);
    }
    fprintf(stdout, "\n");
  }
}

void Print(double& d) { fprintf(stdout, "%.3f", d); }

void extractColumn(Matrix& x, int col, Vector* v) {
  (*v).Dimension(x.rows);
  for (int i = 0; i < x.rows; ++i) {
    (*v)[i] = x[i][col];
  }
}

void extractColumn(Matrix& x, int col, Matrix* m) {
  (*m).Dimension(x.rows, 1);
  for (int i = 0; i < x.rows; ++i) {
    (*m)[i][0] = x[i][col];
  }
}

#endif /* _MATRIXIO_H_ */
