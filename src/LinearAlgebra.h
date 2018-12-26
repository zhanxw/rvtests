#ifndef _LINEARALGEBRA_H_
#define _LINEARALGEBRA_H_

void permute(Vector* v);
void permute(Vector* vec1, Vector* vec2);  // permute both vector.
void centerVector(Vector* v);

inline void permute(Vector* vec) {
  Vector& v = *vec;
  int n = v.Length();
  double tmp;
  for (int i = n - 1; i >= 1; --i) {
    // pick j from 0 <= j <= i
    int j = rand() % (i + 1);
    if (i != j) {
      tmp = v[i];
      v[i] = v[j];
      v[j] = tmp;
    }
  }
}

inline void permute(Vector* vec1, Vector* vec2) {
  Vector& v1 = *vec1;
  Vector& v2 = *vec2;
  int n = v1.Length();
  double tmp;
  for (int i = n - 1; i >= 1; --i) {
    // pick j from 0 <= j <= i
    int j = rand() % (i + 1);
    if (i != j) {
      tmp = v1[i];
      v1[i] = v1[j];
      v1[j] = tmp;

      tmp = v2[i];
      v2[i] = v2[j];
      v2[j] = tmp;
    }
  }
}

inline void centerVector(Vector* v) {
  double avg = v->Average();
  int n = v->Length();
  for (int i = 0; i < n; ++i) {
    (*v)[i] -= avg;
  };
}

inline void centerMatrix(Matrix* v) {
  Matrix& m = *v;
  if (m.rows == 0) return;

  for (int i = 0; i < m.cols; ++i) {
    double s = 0.;
    for (int j = 0; j < m.rows; ++j) {
      s += m(j, i);
    }
    s /= m.rows;
    for (int j = 0; j < m.rows; ++j) {
      m(j, i) -= s;
    }
  }
}

inline void toVector(Matrix& m, int col, Vector* vec) {
  Vector& v = *vec;
  v.Dimension(m.rows);
  for (int i = 0; i < m.rows; ++i) {
    v[i] = m(i, col);
  }
}

/*
 *@return -1: if error happen
 */
inline int corr(Vector& v1, Vector& v2, double* ret) {
  if (v1.Length() != v2.Length()) {
    return -1;
  }
  double s1 = 0.0;
  double s2 = 0.0;
  double s1_2 = 0.0;
  double s2_2 = 0.0;
  double s12 = 0.0;
  int n = v1.Length();
  for (int i = 0; i < n; ++i) {
    s1 += v1[i];
    s1_2 += v1[i] * v1[i];
    s2 += v2[i];
    s2_2 += v2[i] * v2[i];
    s12 += v1[i] * v2[i];
  };
  // fprintf(stderr, "%d %g %g %g %g %g\n", n, s1, s2, s1_2, s2_2, s12);

  double cov = s12 - s1 * s2 / n;
  double var1 = s1_2 - s1 * s1 / n;
  double var2 = s2_2 - s2 * s2 / n;
  double var = (var1 * var2);

  // fprintf(stderr, "%d %g %g %g %g\n", n, cov, var1, var2, var);
  if (var > 1e-30) {
    *ret = cov / sqrt(var);
  } else {
    *ret = 0.0;
  };
  return 0;
}

/*
 *@return -1: if error happen
 */
inline int corr(Matrix& m1, int col1, Vector& v2, double* ret) {
  Vector v1;
  toVector(m1, col1, &v1);
  return (corr(v1, v2, ret));
}

/*
 *@return -1: if error happen
 */
inline int corr(Matrix& m1, int col1, Matrix& m2, int col2, double* ret) {
  Vector v1;
  toVector(m1, col1, &v1);
  Vector v2;
  toVector(m2, col2, &v2);
  return (corr(v1, v2, ret));
}

/**
 * copy m[,col] to @param v
 * @param: 0-based column index
 * @return 0 if OK
 */
inline int extractColumn(Matrix& m, int col, Vector* v) {
  if (col >= m.cols) {
    v = NULL;
    return -1;
  }
  v->Dimension(m.rows);
  for (int i = 0; i < m.rows; ++i) {
    (*v)[i] = m(i, col);
  };
  return 0;
}

/**
 * calculate variance of @param m and @param idx th column
 */
inline double getVariance(Matrix& m, int idx) {
  if (m.rows <= 1) return 0.;

  double s = 0.;
  for (int i = 0; i < m.rows; ++i) {
    s += m(i, idx);
  }
  double avg = s / m.rows;
  s = 0.;
  for (int i = 0; i < m.rows; ++i) {
    s += (m(i, idx) - avg) * (m(i, idx) - avg);
  }
  s /= m.rows;
  return s;
}

inline double getVariance(Vector& v) {
  if (v.Length() <= 1) return 0.;
  const int n = v.Length();
  double s = 0.;
  for (int i = 0; i < n; ++i) {
    s += v[i];
  }
  double avg = s / n;
  s = 0.;
  for (int i = 0; i < n; ++i) {
    s += (v[i] - avg) * (v[i] - avg);
  }
  s /= n;
  return s;
}

inline double getRowVariance(const Matrix& m, int idx) {
  // todo: may need to make it faster, or bypassing this function
  if (m.cols <= 1) return 0.;

  DECLARE_EIGEN_CONST_MATRIX(m, m_e);

  double avg = m_e.row(idx).sum() / m_e.cols();
  return (m_e.row(idx).array() - avg).square().sum() / m_e.cols();

  // double s = 0.;
  // for (int i = 0; i < m.rows; ++i) {
  //   s += m(i,idx);
  // }
  // double avg = s / m.rows;
  // s = 0.;
  // for (int i = 0; i < m.rows; ++i) {
  //   s += (m(i,idx) - avg) * (m(i,idx) - avg);
  // }
  // s /= m.rows;
  // return s;
}

#endif /* _LINEARALGEBRA_H_ */
