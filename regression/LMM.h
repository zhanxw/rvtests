#ifndef LMM_H
#define LMM_H

#include <vector>

class EigenMatrix;
class Vector;
class Matrix;

class LMM {
 public:
  LMM();
  virtual ~LMM();

 public:
  class Impl;
  Impl* impl;

 public:
  // @return 0 when success
  int FitNullModel(Matrix& Xnull, Matrix& y);
  int TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol);
  int AppendKinship(const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  int AppendKinship(const EigenMatrix& kinship);
  int AppendKinship(Matrix& kinship);
  // for Score Test
  double GetUStat() const;
  double GetVStat() const;
  const std::vector<double>& GetSigma2() const;

 public:
  enum Model { EXACT = 0, APPROXIMATE = 1 };
};

#endif /* LMM_H */
