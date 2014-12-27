#include "MultivariateVT.h"
#include "libMvtnorm/mvtnorm.h"

#include <cfloat>
#include <limits>
#include <cmath>
#include <vector>
#include <set>

#include "MatrixOperation.h" // for dump vars

// #define DEBUG
#undef DEBUG

int MultivariateVT::compute(Vector& freq,
                            Matrix& U,
                            Matrix& V) {
  if (freq.Length() != U.rows ||
      freq.Length() != V.rows ||
      U.cols != 1 ||
      V.rows != V.cols) {
    return -1;
  }

  const int numFreq = freq.Length();
  int numKeep = 0;
  std::set<int> skip;
  std::set<int> freqTable; // to avoid float numeric comparison, store maf * 1e6
  maf.Dimension(numFreq);
  for (int i = 0; i < numFreq; ++i) {
    maf[i] = freq[i] < 0.5 ? freq[i] : 1.0 - freq[i];
    if (maf[i] < 1e-10) {
      skip.insert(i);
      continue;
    }
    if (V[i][i] < 1e-10) {
      skip.insert(i);
      continue;
    }
    int mafInt = ceil(maf[i] * 1000000);
    if (freqTable.count(mafInt)) {
      continue;
    }
    freqTable.insert(mafInt);
    numKeep ++;
  }

#ifdef DEBUG
  fprintf(stderr, "numFreq = %d\n", numFreq);  
  fprintf(stderr, "numKeep = %d\n", numKeep);
#endif
  
  if (numKeep == 0) {
    return -1;
  }

  cutoff.Dimension(numKeep);
  // phi[i][j] = 1 means at maf[i] < cutoff[j]
  phi.Dimension(numFreq, numKeep);

  int j = 0;
  for (std::set<int>::const_iterator it = freqTable.begin();
       it != freqTable.end();
       ++it) {
    cutoff[j++] = 1.0 * (*it) / 1000000;
  }

  for (int i = 0; i < numFreq; ++i) {
    for (int j = 0; j < numKeep; ++j) {
      if (skip.count(i) > 0) {
        phi[i][j] = 0.;
        continue;
      }
      if (maf[i] <= cutoff[j]) {
        phi[i][j] = 1.;
      } else {
        phi[i][j] = 0. ;
      }
    }
  }

  // u_phi = t(U) * phi
  Matrix tmp;
  tmp.Transpose(U);
  this->u_phi.Product(tmp, phi);

  // v_phi = t(phi) * v * phi
  Matrix phiT;
  phiT.Transpose(phi);
  tmp.Product(phiT, V);
  this->v_phi.Product(tmp, phi);
  
  int maxIdx = -1;
  double maxVal = -DBL_MAX;
  for (int i = 0; i < numKeep; ++i) {
    double t = fabs(u_phi[0][i] / sqrt(v_phi[i][i]));
    if (t > maxVal) {
      maxIdx = i;
      maxVal = t;
    }
  }

  dumpToFile(freq, "freq");
  dumpToFile(cutoff, "cutoff");    
  dumpToFile(U, "U");
  dumpToFile(V, "V");
  dumpToFile(phi, "phi");
  dumpToFile(u_phi, "u.phi");
  dumpToFile(v_phi, "v.phi");

  this->minMAF = maf.Min();
  this->maxMAF = maf.Max();
  this->optimalMAF = cutoff[maxIdx];
  this->optimalU = u_phi[0][maxIdx];
  this->optimalV = v_phi[maxIdx][maxIdx];
  this->stat = maxVal;
  
  this->optimalNumVar = 0;
  for (int i = 0 ; i < numFreq; ++i) {
    if (phi[i][maxIdx] > 0)
      ++ this->optimalNumVar;
  }
  
  if (this->mvn.getLowerFromCov(maxVal, this->v_phi, &this->pvalue)) {
    this->pvalue = -1; //failed
    return -1;
  }
  this->pvalue = 1.0 - this->pvalue;
  return 0;
}
