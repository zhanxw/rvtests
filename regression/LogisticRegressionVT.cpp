#include "LogisticRegressionVT.h"
#include "MatrixOperation.h"
#include "LogisticRegression.h"
#include "MultivariateNormalDistribution.h" 
#include <Eigen/Dense> // eigen definiation and LDLT

#undef DEBUG
#ifdef DEBUG
#include <iostream>
#include <fstream>

void dump(const Eigen::MatrixXf& m, const char* fn) {
  std::ofstream out(fn);
  out << m;
  out.close();
};
#endif

//////////////////////////////////////////////////
class LogisticRegressionVT::LogisticVTImpl{
 public:
  bool FitNullModel(Matrix& Xnull, Vector& y, int nRound = 100);
  bool TestCovariate(Matrix& Xnull, Vector& y, Matrix& Xcol);
  int GetIndexMax() ;           // return index to the maximum t
  Matrix& GetU() ;
  Matrix& GetV() ;
  Matrix& GetT() ;
  Matrix& GetCov() ;       // return cov(U)
  double GetPvalue() const; 
  double GetEffect(int index) const;
 private:
  void copy(Matrix& m, Eigen::MatrixXf* out) {
    Eigen::MatrixXf& o = *out;
    o.resize(m.rows, m.cols);
    for(int i = 0; i < m.rows; ++i) {
      for (int j = 0; j < m.cols; ++j) {
        o(i,j) = m[i][j];
      }
    }
  }
  void copy(const Eigen::MatrixXf& m, Matrix* out) {
    Matrix& o = *out;
    o.Dimension(m.rows(), m.cols());
    for(int i = 0; i < m.rows(); ++i) {
      for (int j = 0; j < m.cols(); ++j) {
        o[i][j] = m(i,j);
      }
    }
  }
  void copy(Vector& m, Eigen::MatrixXf* out) {
    Eigen::MatrixXf& o = *out;
    int n = m.Length();
    o.resize(n, 1);
    for(int i = 0; i < n; ++i) {
      o(i,0) = m[i];
    }
  }
  
  void rep(double v, int times, std::vector<double>* out) {
    out->resize(times);
    std::fill(out->begin(), out->end(), v);
  }
  /**
   * Convert covariance matrix @param in to correlation matrix, and
   * put lower triangle to @param out
   */
  int cov2cor(const Eigen::MatrixXf& in, std::vector<double>* out) {
    int n = in.rows();
    if ( n != in.cols()) return -1;
    if (n == 0) return -1;
    (*out).resize(n * (n-1)/2);
    int k = 0;
    for (int i = 0 ; i < n; ++i) {
      for (int j = 0; j < i; ++j){
        (*out)[k++] = in(i, j ) / sqrt(in(i,i) * in(j, j));
      }
    }
    return 0;
  }
  /**
   * rescale this->upper, this->lower, and set this->cor
   */
  int makeCov(const Eigen::MatrixXf cov) {
    if (cov.rows() != cov.cols()) return -1;
    if (cov.rows() != (int) upper.size()) return -1;
    if (cov.rows() != (int) lower.size()) return -1;

    for (size_t i = 0; i < upper.size(); ++i) {
      upper[i] /= sqrt(cov(i, i));
      lower[i] /= sqrt(cov(i, i)); 
    }

    cov2cor(cov, &this->cor);
    return 0;
  }
  
  std::vector<double> upper;
  std::vector<double> lower;
  std::vector<double> cor;
  
  Eigen::MatrixXf cov;
  Eigen::MatrixXf y;
  Eigen::MatrixXf geno;
  Eigen::MatrixXf res; // y - y_mle = residual
  Eigen::MatrixXf v;  

  // notations follows Lin AJHG
  Eigen::MatrixXf U; // U_k, k vector
  Eigen::MatrixXf V; // V_k, k vector
  Eigen::MatrixXf t;
  Eigen::MatrixXf Vkk; // MVNormal covariance, k by k matrix
  Eigen::MatrixXf Uk; // U_ki: n by k matrix

  // values will return
  Matrix retU;
  Matrix retV;
  Matrix retT;
  Matrix retCov;
  
  int maxIndex;         // such that abs(t[maxIndex]) is the maximum
  double pvalue;        // store p-values
  LogisticRegression null; // fit null model
  MultivariateNormalDistribution mvnorm;
}; // class LogisticRegressionVT::LogisticVTImpl{

//////////////////////////////////////////////////
// class LogisticRegressionVT //////////////////////////////
LogisticRegressionVT::LogisticRegressionVT() {
  this->impl = new LogisticVTImpl;
};
LogisticRegressionVT::~LogisticRegressionVT() {
  delete this->impl;
}

bool LogisticRegressionVT::FitNullModel(Matrix& Xnull, Vector& y, int nRound){
  return this->impl->FitNullModel(Xnull, y, nRound);
}

bool LogisticRegressionVT::TestCovariate(Matrix& Xnull, Vector& y, Matrix& Xcol) {
  return impl->TestCovariate(Xnull, y, Xcol);
};

int LogisticRegressionVT::GetIndexMax() {
  return this->impl->GetIndexMax();
};

Matrix& LogisticRegressionVT::GetU() {
  return this->impl->GetU();
};

Matrix& LogisticRegressionVT::GetV() {
  return this->impl->GetV();
};

Matrix& LogisticRegressionVT::GetT() {
  return this->impl->GetT();
};

Matrix& LogisticRegressionVT::GetCov() {
  return this->impl->GetCov();
};

double LogisticRegressionVT::GetPvalue() const{
  return this->impl->GetPvalue();
};

double LogisticRegressionVT::GetEffect(int index) const{
  return impl->GetEffect(index);
};

//////////////////////////////////////////////////
// class LogisticRegressionVT::LogisticVTImpl //////////////////////////////
//////////////////////////////////////////////////
bool LogisticRegressionVT::LogisticVTImpl::FitNullModel(Matrix& Xnull, Vector& y, int nRound){
  return (this->null.Fit(Xnull, y));
}

bool LogisticRegressionVT::LogisticVTImpl::TestCovariate(Matrix& Xnull, Vector& yVec, Matrix& Xcol) {
  // const int n = Xnull.rows;
  const int d = Xnull.cols;
  const int k = Xcol.cols;
  copy(Xnull, &cov); // Z = n by d = [z_1^T; z_2^T; ...]
  copy(Xcol, &geno); // S = n by k = [S_1^T, S_2^T, ...]
  copy(yVec, &y);
  copy(this->null.GetPredicted(), &res); // n by 1
  v = res.array() * (1. - res.array());
  res = y - res;
  Eigen::MatrixXf vsz (d, k); // \sum_i v_i S_ki Z_i = n by d matrix
  Eigen::MatrixXf tmp;
  
  // calculate U and V
  const Eigen::MatrixXf& S = geno;
  const Eigen::MatrixXf& Z = cov;
  
  // U = (S.transpose() * (res.asDiagonal())).rowsum();
  U = res.transpose() * S; // 1 by k matrix
 // for (int i = 0; i < d; ++i) {
  //   vsz.col(i) = (Z * (v.array() * S.col().array()).matrix().asDiagonal()).rowsum();
  // }
  vsz = cov.transpose() * v.asDiagonal() * S; // vsz: d by k matrix
  // const double zz = (v.array() * (Z.array().square().rowise().sum()).array()).sum();
  Eigen::MatrixXf zz = cov.transpose() * v.asDiagonal() * cov; // zz: d by d matrix 
  
  // V.size(k, 1);
  //   V(i, 1) = (v.array() * (S.col(i).array().square())).sum() - vsz.row(i).transpose() * vsz.row(i) /   zz;
  // }
  V = (v.asDiagonal() * (S.array().square().matrix())).colwise().sum(); // - // 1 by k
  tmp = ((vsz).array() * (zz.ldlt().solve(vsz)).array()).colwise().sum();
  V -= tmp;

  // V = (v.asDiagonal() * (S.array().square().matrix())).colwise().sum() - ((vsz).array() * (zz.ldlt().solve(vsz)).array()).colwise().sum();

  // Uk is n by k matrix
  // Uk.size(n, k);
  // for (int i = 0; i < k; ++i) {
  //   Uk.col(i) = res * (S.col(i) - vsz.col(i).transpose()) * Z.col(i) / zz;
  // }
  Uk = res.asDiagonal() * (S - Z * zz.ldlt().solve(vsz));
  // Vkk.size(k, k);
  // for (int i = 0; i < k; ++i) {
  //   for (int j = 0; j <= 1; ++j) {
  //     Vkk(i, j) = Uk.col(i) .transpose() * Uk.col(j);
  //   }
  //   if (i != j) {
  //     Vkk(j, i) = Vkk(i, j);
  //   } else {
  //     if (Vkk(i,i) == 0.0) {
  //       return false; // variance term should be larger than zero.
  //     }
  //   }
  //  }
  Vkk = Uk.transpose() * Uk;
  
  Eigen::MatrixXf t = U.array() / V.array().sqrt();
  Eigen::RowVectorXf tmp2 = t.row(0).cwiseAbs();
  tmp2.maxCoeff(&maxIndex);

  rep(-tmp2(maxIndex), k, & lower);
  rep(tmp2(maxIndex), k, & upper);
  makeCov(Vkk);
  if (mvnorm.getBandProbFromCor(k,
                                (double*) lower.data(),
                                (double*) upper.data(),
                                (double*) cor.data(),
                                &pvalue)) {
    fprintf(stderr, "Cannot get MVN pvalue.\n");
    return false;
  }

  copy(U, &retU);
  copy(V, &retV);
  copy(t, &retT);
  copy(Vkk, &retCov);
  pvalue = 1.0 - pvalue;
  return true;
};

int LogisticRegressionVT::LogisticVTImpl::GetIndexMax() {
  return this->maxIndex;
};

Matrix& LogisticRegressionVT::LogisticVTImpl::GetU() {
  return this->retU;
};

Matrix& LogisticRegressionVT::LogisticVTImpl::GetV() {
  return this->retV;
};

Matrix& LogisticRegressionVT::LogisticVTImpl::GetT() {
  return this->retT;
};

Matrix& LogisticRegressionVT::LogisticVTImpl::GetCov() {
  return this->retCov;
};

double LogisticRegressionVT::LogisticVTImpl::GetPvalue() const{
  return this->pvalue;
}; 

double LogisticRegressionVT::LogisticVTImpl::GetEffect(int index) const{
  return U(0, index) / V(0, index);
};
