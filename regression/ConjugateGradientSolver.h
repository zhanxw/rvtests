#ifndef _CONJUGATEGRADIENTSOLVER_H_
#define _CONJUGATEGRADIENTSOLVER_H_

#include "third/eigen/Eigen/Core"

// result = A^(-1) * b
// NOTE: A need to be square
int conjugateSolver(const Eigen::MatrixXf& A, const Eigen::MatrixXf& b,
                    double threshold, Eigen::MatrixXf* result);

// Use conjugate gradient to solve inv(V) * b
// where V = X * X'/M + I * delta
int conjugateSolverBolt(const Eigen::MatrixXf& X, const Eigen::MatrixXf& b,
                        double delta, double threshold,
                        Eigen::MatrixXf* result);

#endif /* _CONJUGATEGRADIENTSOLVER_H_ */
