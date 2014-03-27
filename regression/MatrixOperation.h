#ifndef _MATRIXOPERATION_H_
#define _MATRIXOPERATION_H_

#include "MathMatrix.h"

/**
 * this procedure is essentially:
 * m += v1 * t(v2)
 */
inline void MatrixPlusEqualV1andV2T(Matrix& m, Vector& v1, Vector& v2){
  if (m.rows != v1.Length() || m.cols != v2.Length()){
    fprintf(stderr, "Dimension does not match!");
  };
  for (int i = 0; i < m.rows; i++){
    for (int j = 0; j < m.cols; j++){
      m[i][j] += v1[i] * v2[j];
    }
  }
};

/**
 * this procedure is essentially:
 * m += w * v1 * t(v2)
 */
inline void MatrixPlusEqualV1andV2TWithWeight(Matrix& m, Vector& v1, Vector& v2, double w){
  if (m.rows != v1.Length() || m.cols != v2.Length()){
    fprintf(stderr, "Dimension does not match!");
  };
  for (int i = 0; i < m.rows; i++){
    for (int j = 0; j < m.cols; j++){
      m[i][j] += v1[i] * v2[j] * w;
    }
  }
};

/**
 * Convert column vector to column matrix
 * @param mat = @param v 
 */
inline void copy(Vector& v, Matrix* mat){
  Matrix& m = *mat;
  int n = v.Length();
  m.Dimension(n, 1);
  for (int i = 0; i < n; ++i){
    m[i][0] = v[i];
  }
};

/**
 * Convert matrix @param mat to @parma v
 * Note: only first column is used
 */
inline void copy(Matrix& mat, Vector* v){
  if (mat.cols < 1) return;
  int n = mat.rows;
  v->Dimension(n);
  for (int i = 0; i < n; ++i) {
    (*v)[i] = mat[i][0];
  }
};

inline void print(Vector& v) {
  int n = v.Length();
  fprintf(stderr, "len = %d\n", n);
  for (int i = 0; i < n; ++i){
    fprintf(stderr, "[ %d ] = %g\n", i, v[i]);
  }
};
inline void print(Matrix& mat) {
  int m = mat.rows;
  int n = mat.cols;
  fprintf(stderr, "dim = %d x %d\n", m, n);
  for (int i = 0; i < m; ++i){
    for (int j = 0; j < n; ++j){
      fprintf(stderr, "[ %d, %d ] = %g\t", i, j, mat[i][j]);
    }
    fprintf(stderr, "\n");
  }
};

inline void dumpToFile(Matrix& mat, FILE* fp){
  const int m = mat.rows;
  const int n = mat.cols;
  for (int i = 0; i < m; ++i){
    for (int j = 0; j < n; ++j) {
      if (j)
        fprintf(fp, "\t");
      fprintf(fp, "%g", mat[i][j]);
    }
    fprintf(fp, "\n");
  }
}

inline void dumpToFile(Matrix& mat, const char* fn){
  FILE* fp = fopen(fn, "wt");
  dumpToFile(mat, fp);
  fclose(fp);
};

inline void dumpToFile(Vector& v, const char* fn){
  const int n = v.Length();
  FILE* fp = fopen(fn, "wt");
  for (int i = 0; i < n; ++i){
    fprintf(fp, "%g\n", v[i]);
  }
  fclose(fp);
};

#endif /* _MATRIXOPERATION_H_ */
