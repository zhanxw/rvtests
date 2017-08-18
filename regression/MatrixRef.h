#ifndef _MATRIXREFERENCE_H_
#define _MATRIXREFERENCE_H_

#define REF_TO_EIGEN(matRef, varName)                                   \
  Eigen::Map<Eigen::MatrixXf> varName((matRef).memory_, (matRef).nrow_, \
                                      (matRef).ncol_);

struct FloatMatrixRef {
 public:
  FloatMatrixRef(float* memory, int nrow, int ncol)
      : memory_(memory), nrow_(nrow), ncol_(ncol) {}

 public:
  float* memory_;
  int nrow_;
  int ncol_;
};

#endif /* _MATRIXREFERENCE_H_ */
