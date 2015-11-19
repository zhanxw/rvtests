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
   * res residual of the phenotype
   * v variance
   * X covariate
   * G genotype
   * w weight for G
   * type "C": continuous response; "D": binary response
   * @return 0 when success
   */
  int Fit(Vector& res_G,    // residual under NULL -- may change when permuting
          Vector& v_G,      // variance under NULL -- may change when permuting
          Matrix& X_G,      // covariate
          Matrix& G_G,      // genotype
          Vector& w_G,      // weight
          const char* type  // response type
          );

  double GetPvalue() const;
  double GetQ() const;
  double GetRho() const;

 private:
  // don't copy
  SkatO(const SkatO& s);
  SkatO& operator=(const SkatO& s);
};
#endif
