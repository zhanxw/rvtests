#include "EigenMatrixInterface.h"

#include <fstream>

#include "libsrc/MathMatrix.h"
#include "third/eigen/Eigen/Cholesky"
#include "third/eigen/Eigen/LU"

void G_to_Eigen(Matrix& GM, Eigen::MatrixXf* _EigenM) {
  Eigen::MatrixXf& EigenM = *_EigenM;
  EigenM.resize(GM.rows, GM.cols);
  for (int i = 0; i < GM.rows; i++)
    for (int j = 0; j < GM.cols; j++) EigenM(i, j) = GM[i][j];
}

void G_to_Eigen(Vector& GV, Eigen::VectorXf* _EigenV) {
  Eigen::VectorXf& EigenV = *_EigenV;
  EigenV.resize(GV.Length());
  for (int i = 0; i < GV.Length(); i++) EigenV(i) = GV[i];
}

void Eigen_to_G(const Eigen::MatrixXf& EigenM, Matrix* _GM) {
  Matrix& GM = *_GM;
  GM.Dimension(EigenM.rows(), EigenM.cols());
  for (int i = 0; i < GM.rows; i++)
    for (int j = 0; j < GM.cols; j++) GM[i][j] = EigenM(i, j);
}

void Eigen_to_G(const Eigen::VectorXf& EigenV, Vector* _GV) {
  Vector& GV = *_GV;
  GV.Dimension(EigenV.size());
  for (int i = 0; i < EigenV.size(); i++) GV[i] = EigenV(i);
}

void Eigen_Column_to_G(const Eigen::MatrixXf& EigenM, int colIdx, Vector* _GV) {
  if (colIdx < 0 || colIdx >= EigenM.cols()) return;
  Vector& GV = *_GV;
  GV.Dimension(EigenM.rows());
  for (int i = 0; i < EigenM.rows(); i++) GV[i] = EigenM(i, colIdx);
}

void cbind_G_to_Eigen(Matrix& GM1, Matrix& GM2, Eigen::MatrixXf* EigenM) {
  if (GM1.rows != GM2.rows) return;
  if (!EigenM) return;

  Eigen::MatrixXf& m = *EigenM;
  m.resize(GM1.rows, GM1.cols + GM2.cols);
  for (int i = 0; i < GM1.rows; ++i) {
    for (int j = 0; j < GM1.cols; ++j) {
      m(i, j) = GM1[i][j];
    }
    for (int j = 0; j < GM2.cols; ++j) {
      m(i, j + GM1.cols) = GM2[i][j];
    }
  }
}

//////////////////////////////////////////////////
void G_to_Eigen(Matrix& GM, Eigen::MatrixXd* _EigenM) {
  Eigen::MatrixXd& EigenM = *_EigenM;
  EigenM.resize(GM.rows, GM.cols);
  for (int i = 0; i < GM.rows; i++)
    for (int j = 0; j < GM.cols; j++) EigenM(i, j) = GM[i][j];
}
void G_to_Eigen(Vector& GV, Eigen::VectorXd* _EigenV) {
  Eigen::VectorXd& EigenV = *_EigenV;
  EigenV.resize(GV.Length());
  for (int i = 0; i < GV.Length(); i++) EigenV(i) = GV[i];
}

void Eigen_to_G(const Eigen::MatrixXd& EigenM, Matrix* _GM) {
  Matrix& GM = *_GM;
  GM.Dimension(EigenM.rows(), EigenM.cols());
  for (int i = 0; i < GM.rows; i++)
    for (int j = 0; j < GM.cols; j++) GM[i][j] = EigenM(i, j);
}

void Eigen_to_G(const Eigen::VectorXd& EigenV, Vector* _GV) {
  Vector& GV = *_GV;
  GV.Dimension(EigenV.size());
  for (int i = 0; i < EigenV.size(); i++) GV[i] = EigenV(i);
}

void cbind_G_to_Eigen(Matrix& GM1, Matrix& GM2, Eigen::MatrixXd* EigenM) {
  if (GM1.rows != GM2.rows) return;
  if (!EigenM) return;

  Eigen::MatrixXd& m = *EigenM;
  m.resize(GM1.rows, GM1.cols + GM2.cols);
  for (int i = 0; i < GM1.rows; ++i) {
    for (int j = 0; j < GM1.cols; ++j) {
      m(i, j) = GM1[i][j];
    }
    for (int j = 0; j < GM2.cols; ++j) {
      m(i, j + GM1.cols) = GM2[i][j];
    }
  }
}

void CholeskyInverseMatrix(Matrix& in, Matrix* out) {
  if (in.rows != in.cols) return;

  const int n = in.rows;
  Eigen::MatrixXf x, res;
  G_to_Eigen(in, &x);
  res = x.ldlt().solve(Eigen::MatrixXf::Identity(n, n));
  Eigen_to_G(res, out);
}

double safeSum(const Eigen::MatrixXd& m) {
  const int r = m.rows();
  const int c = m.cols();
  double s = 0.;
  for (int i = 0; i < r; ++i) {
    for (int j = 0; j < c; ++j) {
      if (std::isfinite(m(i, j))) {
        s += m(i, j);
      }
    }
  }
  return s;
}

int matrixRank(Matrix& in) {
  Eigen::MatrixXf x;
  G_to_Eigen(in, &x);

  const int rank = x.fullPivLu().rank();
  return rank;
}

void dumpToFile(const Eigen::MatrixXf& mat, const char* fn) {
  std::ofstream out(fn);
  out << mat.rows() << "\t" << mat.cols() << "\n";
  out << mat;
  out.close();
}

void dumpToFile(const Eigen::Map<Eigen::MatrixXf>& mat, const char* fn) {
  std::ofstream out(fn);
  out << mat.rows() << "\t" << mat.cols() << "\n";
  out << mat;
  out.close();
}
