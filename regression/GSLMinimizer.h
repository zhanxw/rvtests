#ifndef _GSLMINIMIZER_H_
#define _GSLMINIMIZER_H_

#include <gsl/gsl_min.h>

class Minimizer{
public:
  Minimizer();
  ~Minimizer();
  typedef double (Func)(double, void*);
  /**
   * requires f(lowerBound) > f(startValue) < f(upperBound)
   * @return 0 when succeed
   */
  int minimize(gsl_function F, double startValue, double lowerBound, double upperBound);
  double getX() const{return this->finalX;};
  double getY() const{return this->finalY;};
private:
  double finalX; // the point where 
  double finalY;
  
  double epsabs;
  double epsrel;
  int maxIter;
  
  // GSL stuffs
  const gsl_min_fminimizer_type *T;
  gsl_min_fminimizer *s;
};

#endif /* _GSLMINIMIZER_H_ */
