#include "regression/ConjugateGradientSolver.h"

#include "third/eigen/Eigen/Eigen"

// TODO: may compare other cg methods for faster speed,
// e.g. BiCGSTAB
int conjugateSolverBolt(const Eigen::MatrixXf& X, const Eigen::MatrixXf& b,
                        double delta, double threshold,
                        Eigen::MatrixXf* result) {
  Eigen::MatrixXf& x = *result;
  const int M = X.cols();

  x = Eigen::MatrixXf::Zero(b.rows(), 1);
  Eigen::MatrixXf r = b;
  Eigen::MatrixXf r_next;
  Eigen::MatrixXf p = r;
  int k = 1;
  const int kmax = X.cols();
  Eigen::MatrixXf Ap;
  while (k <= kmax) {
#ifdef DEBUG
    fprintf(stderr, "k = %d, time = %s\n", k, currentTime().c_str());
#endif
    // store A * p (or name it V * p)
    Ap = (X * (X.transpose() * p)) / M + p * delta;

    const double r_ss = (r.array().square().sum());
    const double a = r_ss / (p.transpose() * Ap)(0, 0);
    x = x + a * p;
    r_next = r - a * (Ap);
    const double r_next_ss = r_next.array().square().sum();
    if (r_next_ss < threshold) {
      break;
    }
    const double beta = r_next_ss / r_ss;
    p = r_next + beta * p;
    r = r_next;
    k = k + 1;
  }
  return 0;
}

// result = A^(-1) * b
int conjugateSolver(const Eigen::MatrixXf& A, const Eigen::MatrixXf& b,
                    double threshold, Eigen::MatrixXf* result) {
  Eigen::MatrixXf& x = *result;
  x = Eigen::MatrixXf::Zero(b.rows(), 1);
  Eigen::MatrixXf r = b;
  Eigen::MatrixXf r_next;
  Eigen::MatrixXf p = r;
  int k = 1;
  const int kmax = A.cols();
  while (k <= kmax) {
    const double r_ss = (r.array().square().sum());
    const double a = r_ss / (p.transpose() * A * p)(0, 0);
    x = x + a * p;
    r_next = r - a * (A * p);
    const double r_next_ss = r_next.array().square().sum();
    if (r_next_ss < threshold) {
      break;
    }
    const double beta = r_next_ss / r_ss;
    p = r_next + beta * p;
    r = r_next;
    k = k + 1;
  }
  return 0;
}

int conjugateEigen(const Eigen::MatrixXf& A, const Eigen::MatrixXf& b,
                   double threshold, Eigen::MatrixXf* result) {
  Eigen::ConjugateGradient<Eigen::MatrixXf, Eigen::Lower | Eigen::Upper> cg;
  cg.compute(A);
  *result = cg.solve(b);

  return 0;
}
