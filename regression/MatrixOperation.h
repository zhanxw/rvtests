#ifndef _MATRIXOPERATION_H_
#define _MATRIXOPERATION_H_

/**
 * this procedure is essentially:
 * m += v1 * t(v2)
 */
static void MatrixPlusEqualV1andV2T(Matrix& m, Vector& v1, Vector& v2){
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
static void MatrixPlusEqualV1andV2TWithWeight(Matrix& m, Vector& v1, Vector& v2, double w){
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
static void copy(Vector& v, Matrix* mat){
  Matrix& m = *mat;
  int n = v.Length();
  m.Dimension(n, 1);
  for (int i = 0; i < n; ++i){
    mat[i][0] = v[i];
  }
};

#endif /* _MATRIXOPERATION_H_ */
