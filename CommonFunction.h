#ifndef _COMMONFUNCTION_H_
#define _COMMONFUNCTION_H_

#include <gsl/gsl_cdf.h>

template<class T>
class OrderFunction {
public:
OrderFunction(T& t): v(t) {};
  bool operator() (int i, int j)  {
    return v[i] < v[j];
  };
  const T& v;
};

/**
 * @param freq: 0.3, 0.2, 0.1, 0.4
 * will return 0-based rank, does not handle tie
 * @param ord:  2, 1, 0, 3
 *
 * algorithm:
 * make a pair like this: (0.3, 0), (0.2, 1), (0.1, 2), (0.4, 3)
 * sort by first element:
 * (0.1, 2), (0.2, 1), (0.3, 0), (0.4, 3)
 * extract second element:
 * 2, 1, 0, 3
 */
void order(std::vector<double>& freq, std::vector<int>* ord) {
  ord->resize(freq.size());
  for (size_t i = 0; i < freq.size(); ++i)
    (*ord)[i] = i;

  OrderFunction< std::vector<double> > func(freq);
  std::sort(ord->begin(), ord->end(), func);
};

void order(std::vector<int>& freq, std::vector<int>* ord) {
  ord->resize(freq.size());
  for (size_t i = 0; i < freq.size(); ++i)
    (*ord)[i] = i;

  OrderFunction< std::vector<int> > func(freq);
  std::sort(ord->begin(), ord->end(), func);
};

/**
 * Calculate rank, using average for ties.
 */
void calculateRank(std::vector<double>& freq, std::vector<int>* res) {
  std::vector<int> ord;
  order(freq, &ord);
  order(ord, res);
  
  // calculate the mean of rank if there are ties
};

double pnorm(double x) {
  return gsl_cdf_gaussian_P(x, 1.0);
};
double qnorm(double x) {
  return gsl_cdf_gaussian_Pinv(x, 1.0);
};

double calculateMean(const std::vector<double>& v ){
  double s = 0.0;
  if (v.empty()) return 0.0;

  for (unsigned int i = 0; i < v.size(); ++i) {
    s += v[i];
  }
  return s / v.size();
}

// this calculate standard deviation by dividing N
// sd ^ 2 = \frac{1}{N}  (v_i - mean of v_i) ^ 2
double calculateSD(const std::vector<double>& v ){
  double s = 0.0;
  if (v.empty()) return 0.0;

  double mean = calculateMean(v);
  for (unsigned int i = 0; i < v.size(); ++i) {
    const double t = v[i] - mean;
    s += t * t;
  }
  return sqrt(s / v.size());
}

// this calculate standard deviation by dividing N
// sd ^ 2 = \frac{1}{N-1}  (v_i - mean of v_i) ^ 2
double calculateSampleSD(const std::vector<double>& v ){
  double s = 0.0;
  if (v.size() <= 1) return 0.0;

  double mean = calculateMean(v);
  for (unsigned int i = 0; i < v.size(); ++i) {
    const double t = v[i] - mean;
    s += t * t;
  }
  return sqrt(s / ( v.size() - 1));
}

void inverseNormalizeLikeMerlin(std::vector<double>* y){
  if (!y || !y->size()) return;
  const int n = y->size();
  std::vector<int> yRank;
  calculateRank(*y, &yRank);

  // change order to 1-based rank
  for (int i = 0; i < n; ++i ) {
    (*y)[i] = qnorm( ( 0.5 + yRank[i]) / n);
    // fprintf(stderr, "%zu - %g - %d \n", i,  (*y)[i], ord[i]);
  }
  /* double m = calculateMean(*y); */
  /* double sd = calculateSD(*y); */
  // fprintf(stderr, "mean = %g, sd = %g\n", m, sd);
}

// inverse normal a vector AND center it.
void inverseNormalizeLikeR(std::vector<double>* y){
  if (!y || !y->size()) return;
  const int n = y->size();
  std::vector<int> ord;
  order(*y, &ord);

  for (int i = 0; i < n; i++)
    (*y)[i] = ord[i];
  order(*y, &ord);

  double a;
  if ( n <= 10) {
    a = 0.375;
  } else {
    a = 0.5;
  }
  for (int i = 0; i < n ; i++)
    (*y)[i] = qnorm( ( 1.0 + ord[i] - a) / ( n + (1 - a) - a));
}

void standardize(std::vector<double>* y) {
  // center
  size_t n = y->size();
  double m = calculateMean(*y);
  double sd = calculateSD(*y);
  for (size_t i = 0; i < n ; i++) {
    (*y)[i] -= m;
    (*y)[i] /= sd;
  }
  // logger->info("Done: centering to 0.0 and scaling to 1.0 finished.");
};


#endif /* _COMMONFUNCTION_H_ */
