#ifndef __SKAT_H__
#define __SKAT_H__

#include <Eigen/Dense>
#include <fstream>

#include "MixtureChiSquare.h"

class Matrix;
class Vector;

class Skat
{
public:
Skat():pValue(-999) {};
  void Reset() {
    this->pValue = -999.0;
  };
  /**
   * _G suffix: Data structure for Goncalo
   * y phenotype
   * y0 nul model fitted phenotype
   * X covariate
   * v variance
   * G genotype
   * w weight for G
   * @return 0 when success
   */
  int Fit(Vector & res_G,   // residual under NULL -- may change when permuting
          Vector& v_G,      // variance under NULL -- may change when permuting
          Matrix& X_G,      // covariance
          Matrix & G_G,     // genotype
          Vector &w_G);     // weight
  
  double GetQFromNewResidual(Vector & res_G);   // e.g. permuted residual under NULL

  double GetPvalue() const {return this->pValue;};

  double GetQ() const {return this->Q;};

  /* void dump() { */
  /*   dumpToFile(K_sqrt, "out.Ksqrt"); */
  /*   dumpToFile(w_sqrt, "out.Wsqrt"); */
  /*   dumpToFile(P0, "out.P0"); */
  /*   dumpToFile(res, "out.res"); */
  /* }; */
  /* template <class T> */
  /* void dumpToFile(T& v, const char* f){ */
  /*   std::ofstream fout(f); */
  /*   fout << v; */
  /*   fout.close(); */
  /* } */
  
private:
  //Eigen::MatrixXf K;        // G * W * G'
  Eigen::MatrixXf K_sqrt;     // W^{0.5} * G' ----> K = K_sqrt' * K_sqrt
  Eigen::VectorXf w_sqrt;     // W^{0.5} * G' ----> K = K_sqrt' * K_sqrt
  Eigen::MatrixXf P0;         // V - VX ( X' V X)^{-1} X V
  Eigen::VectorXf res;        // residual
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> es;
  // bool hasCache;              // hasCache == true: no need to recalculate P0

  int nPeople;
  int nMarker;
  int nCovariate;
  MixtureChiSquare mixChiSq;

  double pValue;
  double Q;
};

#endif
