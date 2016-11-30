#ifndef _EIGENMATRIX_H_
#define _EIGENMATRIX_H_

#include "third/eigen/Eigen/Core"
/**
 * This class is just a wrapper to Eigen::MatrixXf to speed compling up:
 * other classes just need to forward declare class EigenMatrix
 */
class EigenMatrix {
 public:
  Eigen::MatrixXf mat;
};
#endif /* _EIGENMATRIX_H_ */
