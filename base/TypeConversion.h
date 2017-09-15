#ifndef _TYPECONVERSION_H_
#define _TYPECONVERSION_H_

#include <errno.h>
#include <limits.h>
#include <math.h>    // for HUGE_VALH, HUGE_VALL
#include <stdint.h>  // definition for uint8_t
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <sstream>
#include <vector>

////////////////////////////////////////////////////////////////////////
// define HUGE_VALF, HUGE_VALL for Solaris 10
/* HUGE_VALF is a 'float' Infinity.  */
#ifndef HUGE_VALF
#if defined _MSC_VER
/* The Microsoft MSVC 9 compiler chokes on the expression 1.0f / 0.0f.  */
#define HUGE_VALF (1e25f * 1e25f)
#else
#define HUGE_VALF (1.0f / 0.0f)
#endif
#endif

/* HUGE_VAL is a 'double' Infinity.  */
#ifndef HUGE_VAL
#if defined _MSC_VER
/* The Microsoft MSVC 9 compiler chokes on the expression 1.0 / 0.0.  */
#define HUGE_VAL (1e250 * 1e250)
#else
#define HUGE_VAL (1.0 / 0.0)
#endif
#endif

/* HUGE_VALL is a 'long double' Infinity.  */
#ifndef HUGE_VALL
#if defined _MSC_VER
/* The Microsoft MSVC 9 compiler chokes on the expression 1.0L / 0.0L.  */
#define HUGE_VALL (1e250L * 1e250L)
#else
#define HUGE_VALL (1.0L / 0.0L)
#endif
#endif

// convert double/int/byte to string type
// note: uint8_t by default is treated as char
template <class T>
inline std::string toString(T i) {
  std::stringstream ss;
  ss << i;
  return ss.str();
}

template <class T>
inline std::string toString(const std::set<T>& in, const std::string& sep) {
  if (in.empty()) {
    return "";
  }
  std::stringstream ss;
  typename std::set<T>::const_iterator iter = in.begin();
  for (; iter != in.end(); ++iter) {
    ss << *iter;
    ss << sep;
  }
  std::string ret = ss.str();
  ret.resize(ret.size() - sep.size());
  return ret;
}

template <>
inline std::string toString(const std::set<uint8_t>& in,
                            const std::string& sep) {
  std::stringstream ss;
  std::set<uint8_t>::const_iterator iter = in.begin();
  for (; iter != in.end(); ++iter) {
    ss << (int)*iter;
    ss << sep;
  }
  std::string ret = ss.str();
  ret.resize(ret.size() - sep.size());
  return ret;
}

template <class T>
inline std::string toString(const std::vector<T>& in, const std::string& sep) {
  std::stringstream ss;
  for (size_t i = 0; i != in.size(); ++i) {
    if (i) {
      ss << sep;
    }
    ss << in;
  }
  return ss.str();
}

// convert double/float to string type
// we try to mimic the '%g' in printf
template <class T>
inline std::string floatToString(T i) {
  std::stringstream ss;
  ss.precision(6);
  ss << std::noshowpoint << i;
  return ss.str();
}

template <class T>
inline std::string floatToString(std::vector<T>& i) {
  std::stringstream ss;
  if (i.empty()) {
    return ss.str();
  }

  ss.precision(6);
  ss << std::noshowpoint << i[0];
  for (size_t j = 1; j < i.size(); ++j) {
    ss << ", ";
    ss << std::noshowpoint << i[j];
  }
  return ss.str();
}

// convert int to comma-separated string type
// e.g. -123456 => "-123,456"
std::string toStringWithComma(int in);

// convert std::string to integer
// @return true if conversion succeed
bool str2int(const char* input, int* output);

// convert std::string to integer
// @return true if conversion succeed
inline bool str2int(const std::string& input, int* output) {
  return str2int(input.c_str(), output);
}

// convert std::string to double
// @return true if conversion succeed
inline bool str2double(const char* input, double* output) {
  char* endptr;
  double val;

  errno = 0;
  val = strtod(input, &endptr);

  if (errno == ERANGE) {
#ifndef NDEBUG
    fprintf(stderr, "Over/under flow happened: %s\n", input);
    perror("strtod");
#endif
    return false;
  }
  if (errno == EINVAL) {
    // Ignore error here to avoid displaying:
    // "strtod: Invalid argument" (issue #32)
    // Reason: musl has different implementaiton of strtod,
    // musl set errno = 22 in strtod("NA")
    // glibc set errno = 0 in strtod("NA")
    return false;
  }
  if (errno != 0 && val == 0.) {
#ifndef NDEBUG
    fprintf(stderr, "Unknown conversion error happened: %s\n", input);
    perror("strtod");
#endif
    return false;
  }
  if (endptr == input) {
    // no digits found
    return false;
  }
  *output = val;
  return true;
}
inline bool str2double(const std::string& input, double* output) {
  return str2double(input.c_str(), output);
}

inline int atoi(const std::string& s) {
  int result;
  bool ret = str2int(s.c_str(), &result);
  if (!ret) {
    return 0;
  }
  return result;
};

inline double atof(const std::string& s) {
  double result;
  bool ret = str2double(s.c_str(), &result);
  if (!ret) {
    return 0.0;
  }
  return result;
};

// convert std::string to double
// @return true if conversion succeed
inline bool isdigit(const std::string& s) {
  for (size_t i = 0; i < s.size(); ++i) {
    if (!isdigit(s[i])) return false;
  }
  return true;
}

/**
 * convert @param chrom to integer for easier comparisons
 * leading "chr" will not be considered
 * chr6 -> 6
 * chr6_abc -> 6 + 100
 * chrX -> 23
 * chrY -> 24
 * chrMT -> 25
 * chrOther -> 1000 + ASCII('O')
 * chr -> -1
 */
extern int chrom2int(const std::string& chrom);
#endif /* _TYPECONVERSION_H_ */
