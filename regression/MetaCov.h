#ifndef _METACOV_H_
#define _METACOV_H_

#include <vector>

class EigenMatrix;
class Matrix;
class Vector;

/**
 * This class is to help calculate covariance in family strucutre
 */
class MetaCov{
public: // Make this Impl public to make optimization function easy to write
  class Impl;
  Impl* impl;
public:
  MetaCov();
  ~MetaCov();

  // @return 0 when success
  int FitNullModel(Matrix& Xnull, Matrix& y,
                   const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  // U * ( x - center(x) )
  int TransformCentered(std::vector<double>* x,
                        const EigenMatrix& kinshipU, const EigenMatrix& kinshipS);
  // @param out = sigma2_g * (lambda + delta) ;
  int GetWeight(Vector* out);
};


#endif /* _METACOV_H_ */
