#include <iostream>
#include "regression/MatrixRef.h"
#include "third/eigen/Eigen/Core"

int main() {
  float f[] = {1, 2, 3, 4};
  FloatMatrixRef m(f, 2, 2);
  REF_TO_EIGEN(m, mE);

  assert(mE.rows() == 2);
  assert(mE.cols() == 2);
  assert(mE.sum() == 10);

  std::cout << mE << "\n";

  mE(0, 1) = 6;
  mE(1, 1) = 5;
  assert(f[0] == 1);
  assert(f[1] == 2);
  assert(f[2] == 6);  // default: use column major
  assert(f[3] == 5);

  return 0;
}
