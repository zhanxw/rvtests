#ifndef _VCFVALUE_H_
#define _VCFVALUE_H_

#include <cassert>
#include <string>
#include "VCFBuffer.h"
#include "VCFConstant.h"
#include "base/Exception.h"
#include "base/TypeConversion.h"
#include "base/Utils.h"

class FileWriter;

/**
 * the versatile format to store value for VCF file.
 */
class VCFValue {
 public:
  char* line;  // NOTE: line content may be changed after parsing
  int beg;     // inclusive
  int end;     // exclusive, and beg <= end
  VCFValue() : line(defaultValue), beg(0), end(0) {}
  VCFValue(char* l, int b, int e) : line(l), beg(b), end(e) {}
  static char defaultValue[1];

 public:
  int toInt() const {
    if (!line) return 0;
    return atoi(line + beg);
  }
  void toInt(int* i) const {
    if (!line) {
      *i = 0;
      return;
    }
    *i = atoi(line + beg);
  }
  double toDouble() const {
    if (!line) return 0.0;
    return atof(line + beg);
  }
  void toDouble(double* d) const {
    if (!line) {
      *d = 0.0;
      return;
    }
    *d = strtod(line + beg, 0);
  }
  const char* toStr() const {
    if (!line) return "";
    return (line + beg);
  }
  void toStr(std::string* s) const {
    if (!line) {
      s->resize(0);
      return;
    }
    s->resize(0);
    for (int i = beg; i < end; i++) {
      s->push_back(line[i]);
    }
  }
  void output(FILE* fp) const {
    if (!line) return;
    for (int i = beg; i < end; ++i) {
      fputc(line[i], fp);
    }
  }
  void output(FileWriter* fp) const;

 public:
  // return converted genotypes
  // for multi-allelic genotype such as 0/2, missing will be reported
  int getGenotype() const {
    int g = 0;
    int p = beg;
    if (line[p] == '.') return MISSING_GENOTYPE;
    if (line[p] < '0') {
      REPORT("Wrong genotype detected. [1]");
      return MISSING_GENOTYPE;
    } else {
      g += line[p] - '0';
    }

    if (g > 1) {
      // this is a multi-allelic site
      return MISSING_GENOTYPE;
    }
    p++;
    if (p == end) return g;
    if (line[p] != '|' && line[p] != '/') {
      return MISSING_GENOTYPE;
    }

    p++;
    if (p == end) {
      REPORT("Wrong genotype length = 2");
      return MISSING_GENOTYPE;
    }
    if (line[p] == '.') return MISSING_GENOTYPE;
    if (line[p] < '0') {
      REPORT("Wrong genotype detected. [2]");
    } else {
      const int a2 = line[p] - '0';
      if (a2 > 1) {
        // this is a multi-allelic site
        return MISSING_GENOTYPE;
      }
      g += a2;
    }

    p++;
    if (p != end) {
      return MISSING_GENOTYPE;
    }
    return g;
  }
  /**
   * return male genotype (0 or 2 or missing) that is located in Non-PAR region
   * (1) encoded as haploid, only "0", "1" is allowed and coded as 0, 2
   * (2) encoded as dipolid, only "0|0", "1|1", "0/0", "1/1" is allowed, and
   *     will coded as 0, 2
   * (3) all other genotypes will be coded as missing
   */
  int getMaleNonParGenotype02() const {
    int g = getAllele1();
    if (g == MISSING_GENOTYPE) return MISSING_GENOTYPE;
    if (isHaploid()) {
      if (g == 0) return 0;
      if (g == 1) return 2;
      return MISSING_GENOTYPE;
    }

    // diploid case
    int g2 = getAllele2();
    if (g2 == MISSING_GENOTYPE) return MISSING_GENOTYPE;
    if (g == g2) {
      if (g == 0) return 0;
      if (g == 1) return 2;
    }
    return MISSING_GENOTYPE;
  }
  /**
   * return male genotype (0 or 2 or missing) that is located in Non-PAR region
   * (1) encoded as haploid, only "0", "1" is allowed and coded as 0, 1
   * (2) encoded as dipolid, only "0|0", "1|1", "0/0", "1/1" is allowed, and
   *     will coded as 0, 1
   * (3) all other genotypes will be coded as missing
   */
  int getMaleNonParGenotype01() const {
    int g = getMaleNonParGenotype02();
    if (g == 2) return 1;
    return g;
  }

  /**
   * @return 0 or 1 or 2 or ... as genotype
   */
  int getAllele1() const {
    int g = 0;
    int p = beg;
    if (line[p] == '.') return MISSING_GENOTYPE;
    if (line[p] < '0')
      REPORT("Wrong genotype detected. [1]");
    else
      g += line[p] - '0';
    return g;
  }
  int getAllele2() const {
    int g = 0;
    int p = beg + 2;
    if (p >= end) return MISSING_GENOTYPE;
    if (line[p] == '.') return MISSING_GENOTYPE;
    if (line[p] < '0')
      REPORT("Wrong genotype detected. [2]");
    else
      g += line[p] - '0';
    return g;
  }
  int countAltAllele(int alt) const {
    assert(alt > 0);
    int g = 0;
    int p = beg;
    if (line[p] == '.') return MISSING_GENOTYPE;
    if (line[p] < '0' || line[p] > '9') {
      REPORT("Wrong genotype detected. [2]");
    }
    g += (line[p] - '0' == alt ? 1 : 0);

    p++;
    if (p == end) return g;
    if (line[p] != '|' && line[p] != '/') {
      return MISSING_GENOTYPE;
    }

    p++;
    if (p == end) {
      REPORT("Wrong genotype length = 2");
      return MISSING_GENOTYPE;
    }
    if (line[p] == '.') return MISSING_GENOTYPE;
    if (line[p] < '0' || line[p] > '9') {
      REPORT("Wrong genotype detected. [2]");
    } else {
      g += (line[p] - '0' == alt ? 1 : 0);
    }

    p++;
    if (p != end) {
      return MISSING_GENOTYPE;
    }
    return g;
  }
  int countMaleNonParAltAllele2(int alt) const {
    assert(alt > 0);
    int g = getAllele1();
    if (g == MISSING_GENOTYPE) return MISSING_GENOTYPE;
    assert(0 <= g && g <= 9);
    if (isHaploid()) {
      if (g == alt) {
        return 2;
      } else {
        return 0;
      }
    }

    // diploid case
    int g2 = getAllele2();
    if (g2 == MISSING_GENOTYPE) return MISSING_GENOTYPE;
    if (g == g2) {  // male should have the same alleles
      return (g == alt ? 1 : 0) + (g2 == alt ? 1 : 0);
    }
    return MISSING_GENOTYPE;
  }
  bool isPhased() const {
    if (end - beg != 3) return false;
    int p = beg + 1;
    if (line[p] == '/') return false;
    if (line[p] == '|') return true;
    return false;
  }
  bool isHaploid() const { return (end - beg == 1); }
  bool isMissingGenotype() const {
    if (!line) return true;
    for (int i = beg; i < end; i++) {
      if (line[i] == '.') return true;
    }
    return false;
  }
  /**
   * Just test if the current VCFValue is missing, not if the genotype is
   * missing
   */
  bool isMissing() const { return (beg == end) || line == NULL; }
  /**
   * Use @param s , search from @param beg, till @param c.
   * Store the search results in (*this).
   * @return 0 when success find @param c inside @param s
   * @return 1 when not find c, but reached the end of @param s
   * @return -1: error
   */
  int parseTill(const VCFBuffer& s, const int b, const char c) {
    if (b >= (int)s.size()) return -1;
    this->line = s.getBuffer();
    this->beg = b;
    this->end = b;

    // char* r = strchr(line + beg, c);
    const char* r = ssechr(line + beg, c);
    if (r != NULL) {
      end = r - line;
      return 0;
    } else {
      end = s.size();
      return 1;
    }
  }

  void output(FILE* fp, char c) const {
    for (int i = beg; i <= end; ++i) {
      fputc(line[i], fp);
    }
    fputc(c, fp);
  }
  void output(FileWriter* fp, char c) const;
};

#endif /* _VCFVALUE_H_ */
