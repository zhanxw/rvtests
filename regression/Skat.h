#ifndef __SKAT_H__
#define __SKAT_H__

class Matrix;
class Vector;

class Skat
{
 private:
  class SkatImpl;
  SkatImpl* skatImpl;
 public:
  Skat();
  ~Skat();
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
  int Fit(Vector & res_G,   // residual under NULL -- may change when permuting
          Vector& v_G,      // variance under NULL -- may change when permuting
          Matrix& X_G,      // covariance
          Matrix & G_G,     // genotype
          Vector &w_G);     // weight

  double GetQFromNewResidual(Vector & res_G);   // e.g. permuted residual under NULL

  double GetPvalue() const; //  {return this->pValue;};

  double GetQ() const; // {return this->Q;};

 private:
  // don't copy
  Skat(const Skat& s);
  Skat& operator=(const Skat& s);
};
#endif
