#ifndef _SADDLEPOINTAPPROXIMATION_H_
#define _SADDLEPOINTAPPROXIMATION_H_

#include "third/eigen/Eigen/Dense"

class SaddlePointApproximation {
 public:
  SaddlePointApproximation(const Eigen::MatrixXf& y_,
                           const Eigen::MatrixXf& mu_,
                           const Eigen::MatrixXf& resid_);
  /**
   * @return 0, saddlepoint approximation succeed; -1, failed; -2, approximate
   * works poorly (too far from the mean)
   */

  int calculatePvalue(const Eigen::MatrixXf& g_tilde, float* newPvalue);

 private:
  float CGF_equation(float t, const Eigen::MatrixXf& G,
                     const Eigen::MatrixXf& mu, const Eigen::MatrixXf& y) const;
  float K_function(float t, const Eigen::MatrixXf& G,
                   const Eigen::MatrixXf& mu) const;

  float K_prime_function(float t, const Eigen::MatrixXf& G,
                         const Eigen::MatrixXf& mu) const;

  float K_prime2_function(float t, const Eigen::MatrixXf& G,
                          const Eigen::MatrixXf& mu) const;

 private:
  const Eigen::MatrixXf& y_;
  const Eigen::MatrixXf& mu_;
  const Eigen::MatrixXf& resid_;
};

#endif /* _SADDLEPOINTAPPROXIMATION_H_ */
