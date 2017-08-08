#ifndef _VCFINDIVIDUAL_H_
#define _VCFINDIVIDUAL_H_

#include "VCFBuffer.h"
#include "VCFFunction.h"
#include "VCFValue.h"

#include <vector>

class FileWriter;

// we assume format are always  GT:DP:GQ:GL
class VCFIndividual {
 public:
  // FUNC parseFunction[4];
  VCFIndividual() : fdLen(0) {
    this->include();  // by default, enable everyone
  }
  /**
   * 0-base index for beg and end, e.g.
   *     0 1 2  3
   *     A B C \t
   * beg = 0, end = 3 (line[end] = '\t' or line[end] = '\0')
   *
   * @return
   */
  void parse(const VCFValue& vcfValue) {
    // skip to next
    if (!this->isInUse()) {
      return;
    }

    this->parsed.attach(&vcfValue.line[vcfValue.beg],
                        vcfValue.end - vcfValue.beg);

    size_t idx = 0;
    int beg = 0;
    int ret;
    while (true) {
      if (idx == fdLen) {
        // VCFValue v;
        // fd.push_back(v);
        ++fdLen;
        fd.resize(fdLen);
      }
      VCFValue& v = fd[idx++];
      ret = v.parseTill(this->parsed, beg, ':');
      this->parsed[v.end] = '\0';
      beg = v.end + 1;
      if (ret == 1) {  // last element
        break;
      }
    }
    fdLen = idx;
    if (fdLen == 0) {
      fprintf(stderr, "Empty individual column - very strange!!\n");
      fprintf(stderr, "vcfValue = %s\n", vcfValue.toStr());
    }
  }

  const std::string& getName() const { return this->name; }
  void setName(std::string& s) { this->name = s; }
  void include() { this->inUse = true; }
  void exclude() { this->inUse = false; }
  bool isInUse() { return this->inUse; }

  const VCFValue& operator[](const unsigned int i) const
      __attribute__((deprecated)) {
    if (i >= fdLen) {
      FATAL("index out of bound!");
    }
    return (this->fd[i]);
  }
  VCFValue& operator[](const unsigned int i) __attribute__((deprecated)) {
    if (i >= fdLen) {
      FATAL("index out of bound!");
    }
    return (this->fd[i]);
  }
  /**
   * @param isMissing: index @param i does not exists. Not testing if the value
   * in ith field is missing
   */
  const VCFValue& get(unsigned int i, bool* isMissing) const {
    if (i >= fdLen) {
      *isMissing = true;
      return VCFIndividual::defaultVCFValue;
    }
    *isMissing = this->fd[i].isMissing();
    return (this->fd[i]);
  }
  /**
   * @return VCFValue without checking missingness
   */
  const VCFValue& justGet(unsigned int i) {
    if (i >= fdLen) {
      return VCFIndividual::defaultVCFValue;
    }
    return (this->fd[i]);
  }

  // VCFValue& getSelf() { return this->self; }
  // const VCFValue& getSelf() const { return this->self; }

  size_t size() const { return this->fdLen; }
  /**
   * dump the content of VCFIndividual column
   */
  void output(FILE* fp) const {
    for (size_t i = 0; i < fdLen; ++i) {
      if (i) fputc(':', fp);
      this->fd[i].output(fp);
    }
  }
  void output(FileWriter* fp) const;
  void toStr(std::string* s) const;

 private:
  bool inUse;
  std::string name;  // id name
  // VCFValue self;             // whole field for the individual (unparsed)
  VCFBuffer parsed;          // store parsed string (where \0 added)
  std::vector<VCFValue> fd;  // each field separated by ':', for optimizaiton,
                             // do not use clear(), resize()

  size_t fdLen;  // number of elements in fd
  static VCFValue defaultVCFValue;
  static char defaultValue[2];
};  // end VCFIndividual

#endif /* _VCFINDIVIDUAL_H_ */
