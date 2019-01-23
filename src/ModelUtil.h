#ifndef _MODELUTIL_H_
#define _MODELUTIL_H_

#include <vector>
#include "base/MathMatrix.h"
#include "base/MathVector.h"
#include "base/SimpleMatrix.h"

inline void copy(const std::vector<double>& in, Vector* o) {
  Vector& out = *o;
  out.Dimension(in.size());
  int n = in.size();
  for (int i = 0; i < n; ++i) {
    out[i] = in[i];
  }
}

inline void copyPhenotype(const Matrix& in, Vector* o) {
  Vector& out = *o;
  out.Dimension(in.rows);
  out.data = in.data;
  // for (int i = 0; i < in.rows; ++i) {
  //   out[i] = in(i, 0);
  // }
}

inline void copyPhenotype(const SimpleMatrix& in, Vector* o) {
  Vector& out = *o;
  const int nr = in.nrow();
  out.Dimension(nr);
  for (int i = 0; i < nr; ++i) {
    out[i] = in[i][0];
  }
}

inline void copyVectorToMatrixColumn(Vector& v, SimpleMatrix* out, int column) {
  SimpleMatrix& o = *out;
  const int n = v.Length();
  assert(o.nrow() == n);
  assert(o.ncol() > column);

  for (int i = 0; i < n; ++i) {
    o[i][column] = v[i];
  }
}

/**
 * copy @param in to @param o, with first column being intercept
 */
inline void copyGenotypeWithIntercept(const Matrix& in, Matrix* o) {
  Matrix& out = *o;
  // for (int i = 0; i < in.rows; ++i) {
  //   out(i, 0) = 1.0;
  //   for (int j = 0; j < in.cols; ++j) {
  //     out(i, j + 1) = in(i, j);
  //   }
  // }
  out.Ones(in.rows, 1);
  out.SetColumnLabel(0, "Intercept");

  // for (int i = 0; i < in.cols; ++i)
  //   out.SetColumnLabel(i + 1, in.GetColumnLabel(i));
  out.StackRight(in);
}

/**
 * copy vector of one, @param in and @param cov to @param o (essentially: o =
 * cbind(1, in, cov))
 */
inline void copyGenotypeWithCovariateAndIntercept(const Matrix& in,
                                                  const Matrix& cov,
                                                  Matrix* o) {
  Matrix& out = *o;
  out.Ones(in.rows, 1);
  // for (int i = 0; i < in.rows; ++i) {
  //   out(i, 0) = 1.0;
  // }
  out.SetColumnLabel(0, "Intercept");

  // iter = std::copy(cov.data.begin(), cov.data.end(), iter);
  // for (int j = 0; j < in.cols; ++j) {
  //   for (int i = 0; i < in.rows; ++i) {
  //     out(i, 1 + j) = in(i, j);
  //   }
  //   out.SetColumnLabel(1 + j, in.GetColumnLabel(j));
  // }
  out.StackRight(in);

  // for (int j = 0; j < cov.cols; ++j) {
  //   for (int i = 0; i < cov.rows; ++i) {
  //     out(i, 1 + j + in.cols) = cov(i, j);
  //   }
  //   out.SetColumnLabel(1 + j + in.cols, cov.GetColumnLabel(j));
  // }
  out.StackRight(cov);
}

/**
 * copy intercept and @param cov to @param o with its first column equaling to
 * 1.0
 */
inline void copyCovariateAndIntercept(int n, const Matrix& cov, Matrix* o) {
  if (cov.cols == 0) {
    // (*o).Dimension(n, 1);
    // for (int i = 0; i < n; ++i) {
    //   (*o)(i, 0) = 1.0;
    // }
    (*o).Ones(n, 1);
    (*o).SetColumnLabel(0, "Intercept");
    return;
  }

  // (*o).Dimension(n, 1 + cov.cols);
  // for (int i = 0; i < n; ++i) {
  //   (*o)(i, 0) = 1.0;
  //   for (int j = 0; j < cov.cols; ++j) {
  //     (*o)(i, j + 1) = cov(i, j);
  //   }
  // }
  // (*o).SetColumnLabel(0, "Intercept");
  (*o).Dimension(n, 1);
  std::fill(o->data.begin(), o->data.begin() + cov.rows, 1.0);
  o->SetColumnLabel(0, "Intercept");

  // for (int j = 0; j < cov.cols; ++j) {
  //   (*o).SetColumnLabel(j + 1, cov.GetColumnLabel(j));
  // }
  o->StackRight(cov);
  return;
}

inline void copyCovariateAndIntercept(int n, SimpleMatrix& cov, Matrix* o) {
  const int nr = cov.nrow();
  const int nc = cov.ncol();

  if (nc == 0) {
    (*o).Dimension(n, 1);
    for (int i = 0; i < n; ++i) {
      (*o)(i, 0) = 1.0;
    }

    (*o).SetColumnLabel(0, "Intercept");
    return;
  }
  if (n != nr) {
    assert(false);
    return;
  }

  (*o).Dimension(n, 1 + nc);
  for (int i = 0; i < n; ++i) {
    (*o)(i, 0) = 1.0;
    for (int j = 0; j < nc; ++j) {
      (*o)(i, j + 1) = cov[i][j];
    }
  }
  (*o).SetColumnLabel(0, "Intercept");
  for (int j = 0; j < cov.ncol(); ++j) {
    (*o).SetColumnLabel(j + 1, cov.getColName()[j].c_str());
  }
  return;
}

#endif /* _MODELUTIL_H_ */
