#ifndef _UTILS_H_
#define _UTILS_H_

#include <string.h> // for strlen
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>
#include <string>
#include <cassert>
#include <algorithm>

// transpose matrix with dimension @param nr by @param nc
template<class T>
inline void transposeMatrix(T* m, int nr, int nc) {
  for(int i = 0; i < nr; i++) {
    for (int j = 0; j < i; j++ ){
      std::swap( (*m)[i*nc+j], (*m)[j*nc+i]);
    }
  }
}

// remove the leading 'chr' if any
inline std::string chopChr(const std::string& s) {
  if (s.size() > 3 &&
      (s[0] == 'c' || s[0] == 'C') &&
      (s[1] == 'h' || s[1] == 'H') &&
      (s[2] == 'r' || s[2] == 'R')){
    return s.substr(3);
  }
  return s;
};

// @return true: if @param s has leading "chr", "CHR", "Chr"...
inline bool hasLeadingChr(const std::string& s) {
  if (s.size() > 3 &&
      (s[0] == 'c' || s[0] == 'C') &&
      (s[1] == 'h' || s[1] == 'H') &&
      (s[2] == 'r' || s[2] == 'R')){
    return true;
  }
  return false;
};

// remove the leading and trailing white spaces
inline std::string stringStrip(const std::string& s){
  unsigned int beg = s.find_first_not_of(' ');
  unsigned int end = s.find_last_not_of(' ');
  return s.substr(beg, end-beg);
}
/**
 * split " ",
 *   for "a b", split to "a", "b"
 *   for "a", split to "a"
 *   for "a ", split to "a", ""
 *   for "", split to ""
 */
class StringTokenizer{
 public:
  StringTokenizer(const std::string& i, char token):
      data(i) {
    this->token = token;
    reset();
  }
  StringTokenizer(const std::string& i, const std::string& token):
      data(i) {
    this->token = token;
    reset();
  }
  void reset() {
    this->begin = 0;
    this->end = data.size();
  }
  /**
   * @param piece parsed a piece of string
   * @return true if there are more parsed results
   */
  bool next(std::string* piece) {
    std::string& s = *piece;
    s.clear();
    while (begin <= end) {
      if (begin == end) {
        ++ begin;
        return true;
      }
      
      const char& c = data[begin];
      if (token.find(c) == std::string::npos) {
        // not a token
        s.push_back(c);
        ++begin;
      } else {
        ++begin;
        return begin < end;
      }
    }
    return begin <= end;
  }
 private:
  const std::string& data;
  std::string token;
  size_t begin;
  size_t end;
}; // StringTokenizer

/** tokenize the string
 * @return number of tokens we obtained
 * Special case:
 * For empty input string, we will return 1, and @param result will have only 1 element (the empty string)
 * When delim is empty, we will give warning, return 1, and @param result will have the whole input string
 */
inline int stringTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
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
    if (delim.find(str[i]) != std::string::npos) { // it's a delimeter
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

inline int stringTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
  std::string d(1, delim);
  return (stringTokenize(str, d, result));
};

// pretty much like stringTokenize, but @param result will not contain empty string
inline int stringNaturalTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
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
    if (delim.find(str[i]) != std::string::npos) { // it's a delimeter
      if (s.size()>0){
        result->push_back(s);
        s.clear();
      }
    } else {
      s.push_back(str[i]);
    }
    ++i;
  };
  if (s.size() > 0)
    result->push_back(s);
  return result->size();
};
inline int stringNaturalTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
  std::string d(1, delim);
  return (stringNaturalTokenize(str, d, result));
};

inline void stringJoin(const std::vector<std::string>& array, const char delim, std::string* res){
  res->clear();
  for (size_t i = 0; i < array.size(); i++) {
    if (i) (*res) += delim;
    (*res) += array[i];
  }
};

inline void tolower(std::string* s) {
  for (std::string::iterator i = s->begin();
       i != s->end();
       ++i)
    (*i) = tolower(*i);
};

inline std::string tolower(const std::string& s) {
  std::string ret(s);
  tolower(&ret);
  return ret;
};

inline void toupper(std::string* s) {
  for (std::string::iterator i = s->begin();
       i != s->end();
       ++i)
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


#endif /* _UTILS_H_ */
