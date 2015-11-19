#ifndef __SKATO_H__
#define __SKATO_H__

class Matrix;
class Vector;

class SkatO {
 public:
  class SkatOImpl;
 private:
  SkatOImpl* skatoImpl;

 public:
  SkatO();
  ~SkatO();
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
  int Fit(Vector& res_G,  // residual under NULL -- may change when permuting
          Vector& v_G,    // variance under NULL -- may change when permuting
          Matrix& X_G,    // covariate
          Matrix& G_G,    // genotype
          Vector& w_G);   // weight


  double GetPvalue() const;
  double GetQ() const;
  double GetRho() const; 

 private:
  // don't copy
  SkatO(const SkatO& s);
  SkatO& operator=(const SkatO& s);
};
#endif
