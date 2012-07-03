#ifndef _VCFVALUE_H_
#define _VCFVALUE_H_

#include "Exception.h"
#include "TypeConversion.h"
#include "VCFConstant.h"
#include "VCFBuffer.h"
#include <string>
#include <cassert>

/**
 * the versatile format to store value for VCF file.
 */
class VCFValue{
public:
  const char* line;
  int beg; // inclusive
  int end; // exclusive, and beg <= end
VCFValue():line(NULL), beg(0), end(0){};
public:
  int toInt() const{
    if (!line) return 0;
    return atoi(line+beg);
  };
  void toInt(int* i) const {
    if (!line) {
      *i = 0;
      return;
    }
    *i = atoi(line+beg);
  };
  double toDouble() const {
    if (!line) return 0.0;
    return atof(line+beg);
  };
  void toDouble(double* d) const {
    if (!line) {
      *d = 0.0;
      return ;
    }
    *d = strtod(line+beg, 0);
  };
  const char* toStr() const {
    if (!line) return "";
    return (line + beg);
    /* this->retStr.clear(); */
    /* for (int i = beg; i < end; i++){ */
    /*     this->retStr.push_back(line[i]); */
    /* } */
    /* return (this->retStr.c_str()); */
  };
  void toStr(std::string* s) const {
    if (!line) {
      s->clear();
      return;
    }
    s->clear();
    for (int i = beg; i < end; i++){
      s->push_back(line[i]);
    }
  };
  void output(FILE* fp) const {
    if (!line) return;
    for (int i = beg; i < end; ++i) {
      fputc(line[i], fp);
    }
  };
public:
  // try to convert to genotype
  int getGenotype() const{
    int g = 0;
    int p = beg;
    if (line[p] == '.')
      return MISSING_GENOTYPE;
    if (line[p] < '0')
      REPORT("Wrong genotype detected. [1]");
    else
      g += line[p] - '0';

    p ++ ;
    if (p == end)
      return g;
    if (line[p] != '|' && line[p] != '/') {
      return MISSING_GENOTYPE;
    }

    p ++;
    if (p == end)
      REPORT("Wrong genotype length = 2");
    if (line[p] == '.')
      return MISSING_GENOTYPE;
    if (line[p] < '0')
      REPORT("Wrong genotype detected. [2]");
    else
      g += line[p] - '0';

    p ++;
    if (p != end) {
      return MISSING_GENOTYPE;
    }
    return g;
  };
  /**
   * @return 0 or 1 or 2 as genotype
   */
  int getAllele1() const{
    int g = 0;
    int p = beg;
    if (line[p] == '.')
      return MISSING_GENOTYPE;
    if (line[p] < '0')
      REPORT("Wrong genotype detected. [1]");
    else
      g += line[p] - '0';
    return g;
  };
  int getAllele2() const{
    int g = 0;
    int p = beg + 2;
    if (p >= end)
      return MISSING_GENOTYPE;
    if (line[p] == '.')
      return MISSING_GENOTYPE;
    if (line[p] < '0')
      REPORT("Wrong genotype detected. [2]");
    else
      g += line[p] - '0';
    return g;
  };
  bool isPhased() const{
    if (end - beg != 3) return false;
    int p = beg + 1;
    if (line[p] == '/') return false;
    if (line[p] == '|') return true;
    return false;
  };
  bool isHaploid() const{
    return (end - beg == 1);
  };
  bool isMissingGenotype() const{
    if (!line) return true;
    for (int i = beg; i < end; i++){
      if (line[i] == '.') return true;
    }
    return false;
  };
  /**
   * Just test if the current VCFValue is missing, not if the genotype is missing
   */
  bool isMissing() const{
    return (beg == end) || line == NULL;
  };
  /**
   * Use @param s , search from @param beg, till @param c.
   * Store the search results in (*this).
   * @return 0 when success find @param c inside @param s
   * @return 1 when not find c, but reached the end of @param s
   * @return -1: error
   */
  int parseTill(const VCFBuffer& s, const int b, const char c) {
    if (b >= (int)s.size()) return -1;
    this->line = s.c_str();
    this->beg = b;
    this->end = b;
    while (this->end < (int)s.size()) {
      if (s[this->end] == c) {
        return 0;
      }
      ++this->end;
    }
    return 1;
  };
};

#endif /* _VCFVALUE_H_ */
