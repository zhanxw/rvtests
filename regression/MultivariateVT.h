#ifndef _MULTIVARIATEVT_H_
#define _MULTIVARIATEVT_H_

#include "MathVector.h"
#include "MathMatrix.h"

#include "MultivariateNormalDistribution.h"

// Implement Liu's meta-analysis method
class MultivariateVT{
 public:
  int compute(Vector& freq,
              Matrix& U,
              Matrix& V);
  
  double getOptimalFreq() const {
    return this->optimalFreq;
  }
  double getStat() const {
    return this->stat;
  }
  double getPvalue() const {
    return this->pvalue;
  }
  
 private:
  Vector maf;
  Matrix u;
  Matrix v;
  Matrix phi;
  
  double optimalFreq;
  double stat;
  double pvalue;

  MultivariateNormal mvn;
};

#endif
