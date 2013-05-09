#ifndef _MIXTURECHISQUARE_H_
#define _MIXTURECHISQUARE_H_

#include <cstddef>

class MixtureChiSquare {
public:
MixtureChiSquare():
  sigma(0.0), lim(10000), acc(0.0001) {
    lambda = new double[10];
    noncen = new double[10];
    df = new int[10];
    lambda_size = 0;
    lambda_cap = 10;
  }
  ~MixtureChiSquare(){
    if (lambda_size) {
      delete [] lambda;
      delete [] noncen;
      delete [] df;
      lambda = NULL;
      noncen = NULL;
      df = NULL;
      lambda_size = 0;
      lambda_cap = 0;
    }
  }
  void reset() {
    this->lambda_size = 0;
  };
  void addLambda(double l){
    if (lambda_size +1 == lambda_cap) {
      resize();
    };
    lambda[lambda_size] = l;
    noncen[lambda_size] = 0.0;
    df[lambda_size] = 1;
    ++lambda_size;
  };
  void resize() {
    int newCap = lambda_cap * 2;
    double* newLambda = new double[newCap];
    double* newNonCen = new double[newCap];
    int* newDf = new int[newCap];
    for (int i = 0; i < lambda_cap; ++i){
      newLambda[i] = lambda[i];
      newNonCen[i] = noncen[i];
      newDf[i] = df[i];
    }
    delete[] lambda;
    delete[] noncen;
    delete[] df;
    lambda = newLambda;
    noncen = newNonCen;
    df = newDf;
    lambda_cap = newCap;
  };
  double getPvalue(double Q);
  void dumpLambda() const;
private:
  const double sigma;
  const int lim;
  const double acc;

  // fit in parameters to qf()
  double *lambda;
  double *noncen;
  int *df;
  int lambda_size;
  int lambda_cap;
};

#endif /* _MIXTURECHISQUARE_H_ */
