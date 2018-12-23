#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#include "SaddlePointApproximation.h"

#include "regression/EigenMatrixInterface.h"
#include "regression/MatrixOperation.h"     // dumpToFile
#include "third/gsl/include/gsl/gsl_cdf.h"  // use gsl_cdf_chisq_Q

SaddlePointApproximation::SaddlePointApproximation(const Eigen::MatrixXf& y,
                                                   const Eigen::MatrixXf& mu,
                                                   const Eigen::MatrixXf& resid)
    : y_(y), mu_(mu), resid_(resid) {}

int SaddlePointApproximation::calculatePvalue(const Eigen::MatrixXf& g_tilde,
                                              float* newPvalue) {
  // find root
  // 1. find boundary of the root
  //    first try [0, 5] or [0, -5]
  const float s = (g_tilde.array() * resid_.array()).sum();
  const float var_s =
      (g_tilde.array().square() * (mu_.array() * (1.0 - mu_.array())).abs())
          .sum();
  // if (fabs(s) < 2.0 * sqrt(var_s) &&
  // !std::getenv("BOLTLMM_FORCE_SADDLEPOINT")) {
  //   fprintf(stderr, "Skip saddle point approximation (far from mean, |%g| >
  //   2.0 * sqrt(%g) )!\n", s, var_s);
  //   return -2;
  // }
  float t_grid[2];
  float y_grid[2];
  if (s > 0) {  // \hat{t} > 0
    t_grid[0] = -0.01;
    t_grid[1] = std::min(s / K_prime2_function(0, g_tilde, mu_), (float)5.);

  } else {
    t_grid[0] = 0.01;
    t_grid[1] = std::max(s / K_prime2_function(0, g_tilde, mu_), (float)-5.);
  }
  y_grid[0] = CGF_equation(t_grid[0], g_tilde, mu_, y_);
  y_grid[1] = CGF_equation(t_grid[1], g_tilde, mu_, y_);

  int iter = 0;
  // extend boundary
  // e.g. [0, 5] => [5, 15] => [15, 60]...
  while (y_grid[0] * y_grid[1] > 0) {
    t_grid[0] = t_grid[1];
    y_grid[0] = y_grid[1];
    t_grid[1] *= 4;
    y_grid[1] = CGF_equation(t_grid[1], g_tilde, mu_, y_);

    ++iter;
    fprintf(stderr, "iter %d, %g -> %g \n", iter, t_grid[1], y_grid[1]);

    if (iter > 10) {
      fprintf(stderr,
              "after 10 iteration, still cannot find boundary conditions\n");
      dumpToFile(y_, "tmp.y");
      dumpToFile(g_tilde, "tmp.g_tilde");
      dumpToFile(mu_, "tmp.mu");
      dumpToFile(resid_, "tmp.resid");
      return 1;  // exit(1);
    }
  }
  if (t_grid[0] > t_grid[1]) {
    std::swap(t_grid[0], t_grid[1]);
    std::swap(y_grid[0], y_grid[1]);
  }
  assert(y_grid[0] < 0 && y_grid[1] > 0);  // TODO: make this more robust
  if (y_grid[0] * y_grid[1] > 0) {
    fprintf(stderr, "Wrong boundary conditions!\n");
    // dumpToFile(covEffect, "tmp.covEffect");
    // dumpToFile(xbeta, "tmp.xbeta");
    dumpToFile(g_tilde, "tmp.g_tilde");
    dumpToFile(mu_, "tmp.mu");
    dumpToFile(resid_, "tmp.resid");
    // exit(1);
    return 1;
  }
  float t_new = t_grid[0];
  float y_new;
  iter = 0;
  // stop condition:
  // 1. secant method narrows to a small enough region
  // 2. new propose point has the differnt sign of statistic s. Usually they
  // have the same sign. Unless numerical issue arises.
  while (fabs(t_grid[1] - t_grid[0]) > 1e-3 || t_new * s < 0) {
    t_new = t_grid[0] +
            (t_grid[1] - t_grid[0]) * (-y_grid[0]) / (y_grid[1] - y_grid[0]);

    // switch to bisect?
    const float dist = t_grid[1] - t_grid[0];
    if (t_grid[1] - t_new < 0.1 * dist) {
      t_new = t_grid[0] + (t_grid[1] - t_grid[0]) * 0.8;
    }
    if (t_new - t_grid[0] < 0.1 * dist) {
      t_new = t_grid[0] + (t_grid[1] - t_grid[0]) * 0.2;
    }

    y_new = CGF_equation(t_new, g_tilde, mu_, y_);
    if (y_new == 0) {
      break;
    } else if (y_new > 0) {
      t_grid[1] = t_new;
      y_grid[1] = y_new;
    } else if (y_new < 0) {
      t_grid[0] = t_new;
      y_grid[0] = y_new;
    }
    ++iter;
    fprintf(stderr, "%g -> %g, %g -> %g \n", t_grid[0], y_grid[0], t_grid[1],
            y_grid[1]);

    if (iter > 10) {
      fprintf(stderr,
              "after 10 iteration, secant method still cannot find solution\n");
      dumpToFile(y_, "tmp.y");
      dumpToFile(g_tilde, "tmp.g_tilde");
      dumpToFile(mu_, "tmp.mu");
      dumpToFile(resid_, "tmp.resid");
      break;
    }
  }
  if (fabs(t_new) < 1e-4 && !std::getenv("BOLTLMM_FORCE_SADDLEPOINT")) {
    fprintf(stderr, "Skip saddle point approximation (t is too small: %g)\n",
            t_new);
    return -3;
  }

  // calculate new p_value
  const float K = K_function(t_new, g_tilde, mu_);
  const float K_prime2 = K_prime2_function(t_new, g_tilde, mu_);
  float w = sqrt(2 * (t_new * s - K));
  if (t_new < 0) {
    w = -w;
  }
  const float v = t_new * sqrt(K_prime2);
  assert(v / w > 0);

  float stat = w + log(v / w) / w;
  stat = stat * stat;

  if (t_new * s < K) {
    fprintf(stderr, "Wrong sqrt operand!\n");
    fprintf(stderr,
            "K = %g, K_prime2 = %g, s = %g, t_new = %g, w = %g, v = %g, stat "
            "= %g\n",
            K, K_prime2, s, t_new, w, v, stat);

    dumpToFile(g_tilde, "tmp.g_tilde");
    dumpToFile(mu_, "tmp.mu");
    dumpToFile(resid_, "tmp.resid");
    return 1;
    // exit(1);
  }
  if (std::getenv("BOLTLMM_DUMP_SADDLEPOINT")) {
    fprintf(stderr, "%s:%d dump saddlepoint\n", __FILE__, __LINE__);
    dumpToFile(g_tilde, "tmp.g_tilde");
    dumpToFile(mu_, "tmp.mu");
    dumpToFile(resid_, "tmp.resid");
  }
  fprintf(stderr,
          "K = %g, K_prime2 = %g, s = %g, t = %g, w = %g, v = %g, stat = %g\n",
          K, K_prime2, s, t_new, w, v, stat);
  float& pvalue_ = *newPvalue;
  pvalue_ = gsl_cdf_chisq_Q(stat, 1.0);
  return 0;
}

// CGF - S = K'(t) - S
float SaddlePointApproximation::CGF_equation(float t, const Eigen::MatrixXf& G,
                                             const Eigen::MatrixXf& mu,
                                             const Eigen::MatrixXf& y) const {
  float val = (mu.array() * G.array() /
               ((1.0 - mu.array()) * (-G.array() * t).exp() + mu.array()))
                  .sum() -
              (G.array() * (y).array()).sum();

  static int counter = 0;
  counter++;
  fprintf(stderr, "[%d] t = %g, CGF_equation(t) = %g\n", counter, t, val);

  return val;
}

float SaddlePointApproximation::K_function(float t, const Eigen::MatrixXf& G,
                                           const Eigen::MatrixXf& mu) const {
  // NOTE:
  // calculation of (log(1 - mu.array() + mu.array() * (exp(G.array() * t))))
  // can be very inaccurate,
  // when mu*exp(G*t) is very close to mu, 1 - mu + mu*exp(Gt) can cancel out
  // many useful digits.
  // refer: google "when log1p should be used".
  Eigen::ArrayXf ret = (log1p(
      mu.array() * (G.array() * t).unaryExpr<float (*)(const float)>(&expm1f)));
  ret -= t * (G.array() * mu.array());
  return ret.isFinite().select(ret, 0).sum();
}

float SaddlePointApproximation::K_prime_function(
    float t, const Eigen::MatrixXf& G, const Eigen::MatrixXf& mu) const {
  Eigen::MatrixXf tmp = exp(-G.array() * t);

  return (mu.array() * G.array() /
              ((1.0 + mu.array()) * (-G.array() * t).exp() + mu.array()) -
          G.array() * mu.array())
      .sum();
}

float SaddlePointApproximation::K_prime2_function(
    float t, const Eigen::MatrixXf& G, const Eigen::MatrixXf& mu) const {
  Eigen::ArrayXf denom =
      (1.0 - mu.array()) * (-G.array() * t).exp() + mu.array();

  Eigen::ArrayXf tmp = ((1.0 - mu.array()) * mu.array() * G.array() *
                        G.array() * exp(-G.array() * t) / (denom * denom));
  return tmp.isNaN().select(0, tmp).sum();
  // dumpToFile(tmp, "tmp.tmp");
  //   return ((1.0 - mu.array()) * mu.array() * G.array() * G.array() *
  //           exp(-G.array() * t) / (denom * denom))
  //       .sum();
}
