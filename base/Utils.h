#ifndef _UTILS_H_
#define _UTILS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>  // for strlen

#include <algorithm>
#include <cassert>
#include <string>
#include <vector>

#include "base/CommonFunction.h"  // makeset
#include "base/SimpleString.h"
#include "base/TypeConversion.h"  // toString

/**
 * String related utility functions
 */

// remove the leading 'chr' if any
inline std::string chopChr(const std::string& s) {
  if (s.size() > 3 && (s[0] == 'c' || s[0] == 'C') &&
      (s[1] == 'h' || s[1] == 'H') && (s[2] == 'r' || s[2] == 'R')) {
    return s.substr(3);
  }
  return s;
};

// @return true: if @param s has leading "chr", "CHR", "Chr"...
inline bool hasLeadingChr(const std::string& s) {
  if (s.size() > 3 && (s[0] == 'c' || s[0] == 'C') &&
      (s[1] == 'h' || s[1] == 'H') && (s[2] == 'r' || s[2] == 'R')) {
    return true;
  }
  return false;
};

// remove the leading and trailing white spaces
inline std::string stringStrip(const std::string& s) {
  size_t beg = s.find_first_not_of(' ');
  if (beg == std::string::npos) {
    return "";
  }
  size_t end = s.find_last_not_of(' ');
  return s.substr(beg, end - beg);
}

/**
 * split " ",
 *   for "a b", split to "a", "b"
 *   for "a", split to "a"
 *   for "a ", split to "a" (similar to R strsplit())
 *   for " ", split to ""
 *   for "", split to nothing (length = 0)
 */
class StringTokenizer {
 public:
  StringTokenizer(const std::string& i, char token) : data(&i) {
    this->token.resize(127);  // all ascii characters
    setToken(token);
    reset();
  }
  StringTokenizer(const std::string& i, const std::string& token) : data(&i) {
    this->token.resize(127);  // all ascii characters
    setToken(token);
    reset();
  }
  StringTokenizer(char token) : data(NULL) {
    this->token.resize(127);  // all ascii characters
    setToken(token);
  }
  StringTokenizer(const std::string& token) : data(NULL) {
    this->token.resize(127);  // all ascii characters
    setToken(token);
  }
  /**
   * @param piece parsed a piece of string
   * @return true if there are more parsed results
   */
  bool next(std::string* piece) {
    if (current >= end) {
      return false;
    }

    size_t ptr = current;
    const char* d = data->data();
    while (ptr != end) {
      if (inToken(*(d + ptr))) {
        break;
      } else {
        ++ptr;
      }
    }
    // now ptr is either a delim or the end of the data
    piece->assign(d + current, d + ptr);
    // fprintf(stderr, "piece = %s\n", piece->data());
    current = ptr + 1;
    return true;
  }
  int tokenize(const std::string& in, std::vector<std::string>* result) {
    this->data = &in;
    reset();
    return tokenize(result);
  }
  int tokenize(std::vector<std::string>* result) {
    assert(result);
    // scan to determine size of result
    int n = 0;
    const char* d = this->data->data();
    for (size_t i = 0; i != this->end; ++i) {
      if (inToken(*(d + i))) {
        ++n;
      }
    }
    // store results
    result->resize(n);
    size_t i = 0;
    while (this->next(&((*result)[i]))) {
      ++i;
    }
    return n;
  }
  int naturalTokenize(const std::string& in, std::vector<std::string>* result) {
    this->data = &in;
    reset();
    return naturalTokenize(result);
  }
  int naturalTokenize(std::vector<std::string>* result) {
    assert(result);
    result->resize(0);
    // scan to determine size of result
    int n = 0;
    size_t i = 0;
    size_t last = 0;
    for (i = 0; i != this->end; ++i) {
      if (inToken((*data)[i])) {
        if (last != i) {
          ++n;
          last = i + 1;
        }
      }
    }
    if (last != i) {
      ++n;
    }
    // store results
    result->resize(n);
    i = 0;
    while (this->next(&((*result)[i]))) {
      if (!((*result)[i]).empty()) {
        ++i;
      }
    }
    return n;
  }

 private:
  void reset() {
    this->current = 0;
    this->end = this->data->size();
  }
  void cleanToken() {
    size_t n = token.size();
    for (size_t i = 0; i != n; ++i) {
      token[i] = 0;
    }
  }
  void setToken(const std::string& delim) {
    size_t n = delim.size();
    for (size_t i = 0; i != n; ++i) {
      token[delim[i]] = 1;
    }
  }
  void setToken(const char delim) { token[(int)delim] = true; }
  bool inToken(const char c) const { return token[(int)c]; }

 private:
  const std::string* data;
  std::vector<char> token;
  size_t current;
  size_t end;
};  // StringTokenizer

/** tokenize the string
 * @return number of tokens we obtained
 * Special case:
 * For empty input string, we will return 1, and @param result will have only 1
 * element (the empty string)
 * When delim is empty, we will give warning, return -1, and @param result will
 * have the whole input string
 */
inline int stringTokenize(const std::string& str, const std::string& delim,
                          std::vector<std::string>* result) {
  assert(result);
  result->clear();
  if (!delim.size()) {
    fprintf(stderr, "stringTokenize() using an empty delim");
    result->push_back(str);
    return -1;
  }

  std::string s;
  unsigned int l = str.size();
  unsigned int i = 0;
  while (i < l) {
    if (delim.find(str[i]) != std::string::npos) {  // it's a delimeter
      result->push_back(s);
      s.clear();
    } else {
      s.push_back(str[i]);
    }
    ++i;
  };
  result->push_back(s);
  return result->size();
};

inline int stringTokenize(const std::string& str, const char delim,
                          std::vector<std::string>* result) {
  std::string d(1, delim);
  return (stringTokenize(str, d, result));
};

inline std::vector<std::string> stringTokenize(const std::string& str,
                                               const std::string& delim) {
  std::vector<std::string> result;
  stringTokenize(str, delim, &result);
  return result;
};

inline std::vector<std::string> stringTokenize(const std::string& str,
                                               const char delim) {
  std::vector<std::string> result;
  std::string d(1, delim);
  stringTokenize(str, d, &result);
  return result;
};

// pretty much like stringTokenize, but @param result will not contain empty
// string
inline int stringNaturalTokenize(const std::string& str,
                                 const std::string& delim,
                                 std::vector<std::string>* result) {
  assert(result);
  result->resize(0);
  if (!delim.size()) {
    fprintf(stderr, "stringTokenize() using an empty delim");
    result->push_back(str);
    return -1;
  }
  static std::string ss;
  static SimpleString s(4096);
  s.resize(0);
  unsigned int l = str.size();
  unsigned int i = 0;
  while (i < l) {
    if (delim.find(str[i]) != std::string::npos) {  // it's a delimeter
      if (!s.empty()) {
        ss = s.data();
        result->push_back(ss);
        s.resize(0);
      }
    } else {
      s.append(str[i]);
    }
    ++i;
  };
  if (!s.empty()) {
    ss = s.data();
    result->push_back(ss);
  }
  return result->size();
};

inline int stringNaturalTokenize(const std::string& str, const char delim,
                                 std::vector<std::string>* result) {
  std::string d(1, delim);
  return (stringNaturalTokenize(str, d, result));
};

inline int stringTokenize(const std::string& str, const std::string& delim,
                          const std::string& leftQuote,
                          const std::string& rightQuote,
                          std::vector<std::string>* result) {
  // check parameters
  int table[128] = {0};
  for (size_t i = 0; i != delim.size(); ++i) {
    table[(int)delim[i]]++;
  }
  for (size_t i = 0; i != leftQuote.size(); ++i) {
    table[(int)leftQuote[i]]++;
  }
  for (size_t i = 0; i != rightQuote.size(); ++i) {
    table[(int)rightQuote[i]]++;
  }
  for (int i = 0; i < 128; ++i) {
    if (table[i] > 1) {
      fprintf(stderr, "Duplicated character found [ %c ]!\n", (char)i);
      return -1;
    }
  }

  // handle special inputs
  std::string s;
  result->clear();
  if (!delim.size()) {
    result->push_back(s);
    return -1;
  }
  if (!str.size()) {
    result->push_back(s);
    return -1;
  }

  // normal workflow
  size_t len = str.size();
  int state = 0;
  for (size_t i = 0; i != len; ++i) {
    const char c = str[i];
    if (delim.find(c) != std::string::npos) {  // a delim
      if (state == 0) {
        result->push_back(s);
        s.clear();
      } else {
        s.push_back(c);
      }
    } else if (leftQuote.find(c) != std::string::npos) {  // left quote
      state++;
      s.push_back(c);
    } else if (rightQuote.find(c) != std::string::npos) {  // right quote
      state--;
      s.push_back(c);
      if (state < 0) {
        fprintf(stderr, "Unbalanced quotes: [ %s ]\n", str.c_str());
        return -1;
      }
    } else {  //
      s.push_back(c);
    }
  }
  if (state != 0) {
    fprintf(stderr, "Unbalanced quotes: [ %s ]\n", str.c_str());
    return -1;
  }
  if (s.size()) {
    result->push_back(s);
  }
  return result->size();
}

inline void stringJoin(const std::vector<std::string>& array, const char delim,
                       std::string* res) {
  res->clear();
  for (size_t i = 0; i < array.size(); i++) {
    if (i) (*res) += delim;
    (*res) += array[i];
  }
};

inline std::string stringJoin(const std::vector<std::string>& array,
                              const char delim) {
  std::string res;
  stringJoin(array, delim, &res);
  return res;
};

inline void tolower(std::string* s) {
  for (std::string::iterator i = s->begin(); i != s->end(); ++i)
    (*i) = tolower(*i);
};

inline std::string tolower(const std::string& s) {
  std::string ret(s);
  tolower(&ret);
  return ret;
};

inline void toupper(std::string* s) {
  for (std::string::iterator i = s->begin(); i != s->end(); ++i)
    (*i) = toupper(*i);
};

inline std::string toupper(const std::string& s) {
  std::string ret(s);
  toupper(&ret);
  return ret;
};

inline bool isInteger(const double m) {
  double d;
  return (modf(m, &d) == 0.0);
}

/**
 * print out the content for debug only
 */
inline void dumpStringVector(const std::vector<std::string> s) {
  for (unsigned int i = 0; i < s.size(); i++) {
    fprintf(stdout, "%u: %s\n", i, s[i].c_str());
  }
};

/**
 * @return true if @param s ends with @param tail
 */
inline bool endsWith(const std::string& s, const std::string& tail) {
  if (s.size() < tail.size()) return false;
  size_t l = tail.size();
  size_t idx = s.size() - l;
  for (size_t i = 0; i != l; ++i) {
    if (s[idx + i] != tail[i]) {
      return false;
    }
  }
  return true;
}

extern char const* ssechr(char const* s, char ch);

#endif /* _UTILS_H_ */
