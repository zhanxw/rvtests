#ifndef _COMMONFUNCTION_H_
#define _COMMONFUNCTION_H_

#include <algorithm>
#include <map>
#include <set>
#include "base/Indexer.h"
#include "base/Utils.h"
#include "third/gsl/include/gsl/gsl_cdf.h"

extern int stringTokenize(const std::string& str, const char delim,
                          std::vector<std::string>* result);

//////////////////////////////////////////////////
// Statistics functions
//////////////////////////////////////////////////
template <class T>
class OrderFunction {
 public:
  OrderFunction(T& t) : v(t) {}
  bool operator()(int i, int j) { return v[i] < v[j]; }
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
inline void order(std::vector<double>& in, std::vector<int>* ord) {
  ord->resize(in.size());
  for (size_t i = 0; i < in.size(); ++i) (*ord)[i] = i;

  OrderFunction<std::vector<double> > func(in);
  std::sort(ord->begin(), ord->end(), func);
}

inline void order(std::vector<int>& in, std::vector<int>* ord) {
  ord->resize(in.size());
  for (size_t i = 0; i < in.size(); ++i) (*ord)[i] = i;

  OrderFunction<std::vector<int> > func(in);
  std::sort(ord->begin(), ord->end(), func);
}

/**
 * Calculate rank, using average for ties.
 * Rank starts at 0.
 * For tie rank, rank = (rank_min + rank_high) / 2
 */
inline void calculateRank(std::vector<double>& in, std::vector<double>* out) {
  std::map<double, double> counter;
  for (size_t i = 0; i != in.size(); ++i) {
    counter[in[i]]++;
  }
  int low = 0;  // start with lowest rank
  int high = 0;
  std::map<double, double>::iterator it = counter.begin();
  for (; it != counter.end(); ++it) {
    high = low + it->second;
    if (it->second == 1) {  // not a tie
      it->second = low;
    } else {  // encounter tie
      it->second = 0.5 * (low + (low + it->second - 1.0));
    }
    low = high;
  }
  out->resize(in.size());
  for (size_t i = 0; i != in.size(); ++i) {
    (*out)[i] = counter[in[i]];
  }
}

inline double pnorm(double x) { return gsl_cdf_ugaussian_P(x); }
inline double qnorm(double x) { return gsl_cdf_ugaussian_Pinv(x); }
inline double qnorm(double x, double sigma) {
  return gsl_cdf_gaussian_Pinv(x, sigma);
}

inline double calculateMean(const std::vector<double>& v) {
  double s = 0.0;
  if (v.empty()) return 0.0;

  for (unsigned int i = 0; i < v.size(); ++i) {
    s += v[i];
  }
  return s / v.size();
}

inline void centerVector(std::vector<double>* v) {
  if (v->empty()) return;
  double m = calculateMean(*v);
  for (size_t i = 0; i < v->size(); ++i) {
    (*v)[i] -= m;
  }
}

// this calculate standard deviation by dividing N
// sd ^ 2 = \frac{1}{N}  (v_i - mean of v_i) ^ 2
inline double calculateSD(const std::vector<double>& v) {
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
inline double calculateSampleSD(const std::vector<double>& v) {
  double s = 0.0;
  if (v.size() <= 1) return 0.0;

  double mean = calculateMean(v);
  for (unsigned int i = 0; i < v.size(); ++i) {
    const double t = v[i] - mean;
    s += t * t;
  }
  return sqrt(s / (v.size() - 1));
}

inline void inverseNormalizeLikeMerlin(std::vector<double>* y) {
  if (!y || !y->size()) return;
  const int n = y->size();
  std::vector<double> yRank;
  calculateRank(*y, &yRank);

  // change order to 1-based rank
  for (int i = 0; i < n; ++i) {
    (*y)[i] = qnorm((0.5 + yRank[i]) / n);  // rank start with 0, so + 0.5 here
    // fprintf(stderr, "%zu - %g - %d \n", i,  (*y)[i], ord[i]);
  }
  /* double m = calculateMean(*y); */
  /* double sd = calculateSD(*y); */
  // fprintf(stderr, "mean = %g, sd = %g\n", m, sd);
}

// inverse normal a vector AND center it.
inline void inverseNormalizeLikeR(std::vector<double>* y) {
  if (!y || !y->size()) return;
  const int n = y->size();
  std::vector<int> ord;
  order(*y, &ord);

  for (int i = 0; i < n; i++) (*y)[i] = ord[i];
  order(*y, &ord);

  double a;
  if (n <= 10) {
    a = 0.375;
  } else {
    a = 0.5;
  }
  for (int i = 0; i < n; i++)
    (*y)[i] = qnorm((1.0 + ord[i] - a) / (n + (1 - a) - a));
}

inline void zero(std::vector<double>* y) {
  std::fill(y->begin(), y->end(), 0.0);
}

inline void standardize(std::vector<double>* y) {
  if (!y) return;
  const size_t n = y->size();
  // center
  double m = calculateMean(*y);
  double sd = calculateSD(*y);
  if (sd == 0.0) {
    zero(y);
    return;
  }

  for (size_t i = 0; i != n; i++) {
    (*y)[i] -= m;
    (*y)[i] /= sd;
  }
  // logger->info("Done: centering to 0.0 and scaling to 1.0 finished.");
}

/**
 * Look up an element @param elem from @param vec
 * @return index if found, or -1 if not.
 */
inline int which(const std::vector<std::string>& vec, const char* elem) {
  for (size_t i = 0; i != vec.size(); ++i) {
    if (vec[i] == elem) {
      return i;
    }
  }
  return -1;
}

/**
 * Construct a indice @param index for @param a such that a[index] == b
 * NOTE: similar to mathc() in R.
 * @param nomatch, when an element in @param a are not in @param b, this value
 * will be used as an index
 * @return number of elements that are matched
 * e.g. a = ["D", "B", "A"] b = ["A", "B", "C"],
 * match(a, b, index) will return [-1, 1, 0]
 */
inline int match(const std::vector<std::string>& a,
                 const std::vector<std::string>& b, std::vector<int>* index,
                 int nomatch = -1) {
  int numMatched = 0;
  Indexer idx(b);
  index->reserve(a.size());
  for (size_t i = 0; i != a.size(); ++i) {
    index->push_back(idx[a[i]]);
    if (index->back() == -1) {
      index->back() = nomatch;
    } else {
      numMatched++;
    }
  }

  return numMatched;
}

//////////////////////////////////////////////////
// STL based functions related to set
//////////////////////////////////////////////////

/**
 * Convert a @param string separated by @param sep to set (stored in @param s)
 */
inline void makeSet(const std::string& str, char sep,
                    std::set<std::string>* s) {
  s->clear();
  if (str.empty()) return;

  std::vector<std::string> fd;
  stringTokenize(str, sep, &fd);
  for (size_t i = 0; i < fd.size(); i++) s->insert(fd[i]);
}

template <class T>
inline void makeSet(const std::vector<T>& in, std::set<T>* s) {
  s->clear();
  if (in.empty()) return;

  for (size_t i = 0; i < in.size(); i++) s->insert(in[i]);
}

template <class T>
inline std::set<T> makeSet(const std::vector<T>& in) {
  std::set<T> s;
  s.clear();
  for (size_t i = 0; i < in.size(); i++) s.insert(in[i]);
  return s;
}

template <class T>
inline void makeCounter(const std::vector<T>& in, std::map<T, int>* s) {
  s->clear();
  if (in.empty()) return;

  for (size_t i = 0; i < in.size(); i++) {
    (*s)[in[i]]++;
  }
}

/**
 * make a vector to map.
 *   when there is no duplciation: key is vector[i], value is i
 *   when there is duplication: key is vector[i], value is the index of the
 * first appearance of vector[i]
 */
inline void makeMap(const std::vector<std::string>& in,
                    std::map<std::string, int>* s) {
  s->clear();
  if (in.empty()) return;

  for (size_t i = 0; i < in.size(); i++) {
    if (s->find(in[i]) != s->end()) continue;
    (*s)[in[i]] = i;
  }
}

inline void makeMap(const std::string& in, std::map<char, int>* s) {
  s->clear();
  if (in.empty()) return;

  for (size_t i = 0; i < in.size(); i++) {
    if (s->find(in[i]) != s->end()) continue;
    (*s)[in[i]] = i;
  }
}

inline std::map<char, int> makeMap(const std::string& in) {
  std::map<char, int> s;
  makeMap(in, &s);
  return s;
}

/**
 *  Test whether all elements in @param x are unique
 */
inline bool isUnique(const std::vector<std::string>& x) {
  std::set<std::string> s;
  for (size_t i = 0; i < x.size(); i++) {
    s.insert(x[i]);
    if (s.size() != i + 1) {
      return false;
    }
  }
  return true;
}

/**
    * Remove i th element from @param val where i is stored in @param index
    * @return number of elements removed
    *
    * NOTE: unless static function, template functions should not be in .cpp
  *files
    */
template <typename T, typename A>
inline int removeByIndex(const std::vector<int>& index,
                         std::vector<T, A>* val) {
  if (index.empty()) return 0;

  std::set<int> idx(index.begin(), index.end());

  int nRemoved = 0;
  size_t last = 0;
  const size_t n = val->size();
  for (size_t i = 0; i < n; ++i) {
    if (idx.count(i)) {
      ++nRemoved;
      continue;
    }
    if (last != i) {
      (*val)[last] = (*val)[i];
    }
    ++last;
  }
  val->resize(last);
  return nRemoved;
}

template <class T>
struct FuncInSet {
  FuncInSet(const std::set<T>& s) : data(s) {}
  bool operator()(T& s) const {
    if (this->data.count(s)) {
      return true;
    }
    return false;
  }
  const std::set<T>& data;
};
template <class T>
struct FuncNotInSet {
  FuncNotInSet(const std::set<T>& s) : data(s) {}
  bool operator()(T& s) const {
    if (this->data.count(s)) {
      return false;
    }
    return true;
  }
  const std::set<T>& data;
};

template <class T>
inline int intersect(const std::vector<T>& a, std::vector<T>* b) {
  if (!b) return -1;
  std::set<T> s(a.begin(), a.end());
  FuncNotInSet<T> func(s);
  b->erase(std::remove_if(b->begin(), b->end(), func), b->end());
  return 0;
}

template <class T>
inline std::vector<T> intersect(const std::vector<T>& a,
                                const std::vector<T>& b) {
  std::vector<T> ret(b);
  intersect(a, &ret);
  return ret;
}

/**
 * Test all elements in @param a exist in @param b
 */
template <class T>
bool isSubset(const std::vector<T>& a, const std::vector<T>& b) {
  std::set<T> s;
  makeSet(b, &s);
  for (typename std::vector<T>::const_iterator it = a.begin(); it != a.end();
       ++it) {
    if (s.count(*it) == 0) {
      return false;
    }
  }
  return true;
}

/**
 * Remove elements from @param b if they exist in @param a, i.e. set(b) - set(a)
 */
template <class T>
inline int remove(const std::vector<T>& a, std::vector<T>* b) {
  if (!b) return -1;
  std::set<T> s(a.begin(), a.end());
  FuncInSet<T> func(s);
  b->erase(std::remove_if(b->begin(), b->end(), func), b->end());
  return 0;
}

/**
 * Return set(a) - set(b)
 */
template <class T>
inline std::vector<T> setSubtract(const std::vector<T>& a,
                                  const std::vector<T>& b) {
  std::vector<T> s(a);
  remove(b, &s);
  return s;
}

template <class T>
inline int dedup(std::vector<T>* a) {
  if (!a) return -1;
  std::set<T> s;
  size_t idx = 0;
  size_t n = a->size();
  for (size_t i = 0; i != n; ++i) {
    if (s.count((*a)[i])) {
      continue;
    }
    s.insert((*a)[i]);
    if (idx != i) {
      (*a)[idx] = (*a)[i];
    }
    idx++;
  }

  a->resize(idx);
  return 0;
}

template <class T>
inline int extend(const std::vector<T>& vec, std::vector<T>* a) {
  if (!a) return -1;
  size_t n = vec.size();
  for (size_t i = 0; i != n; ++i) {
    a->push_back(vec[i]);
  }
  return 0;
}

#if 0
template <class T>
int which(const std::vector<T>& vec, T& elem) {
  for (size_t i = 0; i != vec.size(); ++i) {
    if (vec[i] == elem) {
      return i;
    }
  }
  return -1;
}
#endif
#endif /* _COMMONFUNCTION_H_ */
