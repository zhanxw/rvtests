#include "MultivariateVT.h"
#include "libMvtnorm/mvtnorm.h"

#include <cfloat>
#include <limits>
#include <cmath>
#include <vector>
#include <map>
// MultivariateVT::MultivariateVT(Vector& freq,
//                                Matrix& U,
//                                Matrix& V):
//     freq(freq),
//     U(U),
//     V(V)
// {
// }

int MultivariateVT::compute(Vector& freq,
                            Matrix& U,
                            Matrix& V) {
  if (freq.Length() != U.rows ||
      freq.Length() != V.rows ||
      U.cols != 1 ||
      V.rows != V.cols) {
    return -1;
  }

  const int n = freq.Length();
  int numKeep = 0;
  std::map<double, std::vector<int> > freqTable;
  for (int i = 0; i < n; ++i) {
    double f  = freq[i] < 0.5 ? freq[i] : 1.0 - freq[i];
    if (f < 1e-10) continue;
    if (V[i][i] < 1e-10) continue;
    freqTable[f].push_back(i);
    numKeep ++;
  }

  maf.Dimension(numKeep);
  u.Dimension(numKeep, 1);
  v.Dimension(numKeep, numKeep);
  phi.Dimension(numKeep, numKeep);

  int i = 0;
  for (std::map<double, std::vector<int> >::const_iterator it = freqTable.begin();
       it != freqTable.end();
       ++it) {
    double f = it->first;
    const std::vector<int>& idx = it->second;

    maf[i] = f;
    for (int j = 0; j < (int) idx.size(); ++j) {
      u[i][0] = U[idx[j]][0];

      for (int k = 0; k < n; ++k ){
        v[k][i] = V[k][idx[j]];
      }
    }
    i ++;
  }
  for (int i = 0; i < numKeep; ++i) {
    for (int j = 0; j < numKeep; ++j) {
      if (maf[i] <= maf[j]) {
        phi[i][j] = 1.0;
      } else {
        phi[i][j] = 0.0;
      }
    }
  }

  int maxIdx = -1;
  double maxVal = -DBL_MAX;
  for (int i = 0; i < numKeep; ++i) {
    double t = fabs(u[i][0] / v[i][i]);
    if (t > maxVal) {
      maxIdx = i;
      maxVal = t;
    }
  }

  Matrix phiT;
  phiT.Transpose(phi);
  Matrix tmp;
  tmp.Product(phi, v);
  Matrix cov;
  cov.Product(tmp, phiT);

  this->optimalFreq = maf[maxIdx];
  this->stat = maxVal;

  if (this->mvn.getUpperFromCov(fabs(this->stat), cov, &this->pvalue)) {
    this->pvalue = -1; //failed
    return -1;
  }

  return 0;
}
