#ifndef _MULTIVARIATEVT_H_
#define _MULTIVARIATEVT_H_

#include "MathVector.h"
#include "MathMatrix.h"

#include "MultivariateNormalDistribution.h"

// Implement Liu's meta-analysis method
class MultivariateVT{
 public:
  /**
   * @param freq
   * @param U n by 1 column matrix
   * @param V n by n square matrixde
   */
  int compute(Vector& freq,
              Matrix& U,
              Matrix& V);
  
  double getMinMAF() const {
    return this->minMAF;
  }
  double getMaxMAF() const {
    return this->maxMAF;
  }
  double getOptimalMAF() const {
    return this->optimalMAF;
  }
  int getOptimalNumVar() const {
    return this->optimalNumVar;
  }
  double getOptimalU() const {
    return this->optimalU;
  }
  double getOptimalV() const {
    return this->optimalV;
  }
  double getStat() const {
    return this->stat;
  }
  double getPvalue() const {
    return this->pvalue;
  }
  
 private:
  Vector maf;
  Vector cutoff;
  Matrix u_phi; // U * phi
  Matrix v_phi; // phi' * V * phi
  Matrix phi;

  double minMAF;  
  double maxMAF;
  double optimalMAF;
  int optimalNumVar;  
  double optimalU;
  double optimalV;
  double stat;
  double pvalue;

  MultivariateNormalDistribution mvn;
};

#endif
