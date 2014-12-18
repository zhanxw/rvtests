#ifndef __FAMSKAT_H__
#define __FAMSKAT_H__

class Matrix;
class Vector;
class EigenMatrix;

class FamSkat
{
 private:
  class FamSkatImpl;
  FamSkatImpl* skatImpl;
 public:
  FamSkat();
  ~FamSkat();
  void Reset();
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
  int FitNullModel(Matrix& Xnull, Matrix& y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  int TestCovariate(Matrix& Xnull, Matrix& y, Matrix& Xcol, Vector& weight,
                    const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);

  // not implemented yet,
  // but doable using parametric bootstrap method
  // double FamSkat::GetQFromPermutation();

  double GetPvalue() const; //  {return this->pValue;};

  double GetQ() const; // {return this->Q;};

 private:
  // don't copy
  FamSkat(const FamSkat& s);
  FamSkat& operator=(const FamSkat& s);
};
#endif
