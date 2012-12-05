#ifndef _VCFINDIVIDUAL_H_
#define _VCFINDIVIDUAL_H_

#include "VCFFunction.h"
#include "VCFValue.h"
#include "VCFBuffer.h"

#include <vector>

// we assume format are always  GT:DP:GQ:GL
class VCFIndividual{
public:

  // FUNC parseFunction[4];
  VCFIndividual(){
    this->include();  // by default, enable everyone
  };
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

    this->self = vcfValue;
    this->parsed = vcfValue.toStr();

    // need to consider missing field
    this->fd.clear();

    VCFValue v;
    int beg = 0;
    int ret;
    while( (ret = v.parseTill(this->parsed, beg, ':')) == 0){
      this->parsed[v.end] = '\0';
      beg = v.end + 1;
      fd.push_back(v);
    }
    if (ret == 1) {
      this->parsed[v.end] = '\0';
      fd.push_back(v);
    }
    if (fd.size() == 0) {
      fprintf(stderr, "Empty individual column - very strange!!\n");
      fprintf(stderr, "vcfValue = %s\n", vcfValue.toStr());
    }
  };

  const std::string& getName() const {return this->name;};
  void setName(std::string& s) {this->name = s;};
  void include() {this->inUse = true;};
  void exclude() {this->inUse = false;};
  bool isInUse() {return this->inUse;};

  const VCFValue& operator [] (const unsigned int i) const __attribute__ ((deprecated)) {
    if (i >= fd.size()){
      FATAL("index out of bound!");
    }
    return (this->fd[i]);
  };
  VCFValue& operator [] (const unsigned int i) __attribute__ ((deprecated)) {
    if (i >= fd.size()){
      FATAL("index out of bound!");
    }
    return (this->fd[i]);
  };
  const VCFValue& get(unsigned int i, bool* isMissing) const {
    if (i >= fd.size()) {
      *isMissing = true;
      return VCFIndividual::defaultVCFValue;
    }
    *isMissing = this->fd[i].isMissing();
    return (this->fd[i]);
  };
  /**
   * @return VCFValue without checking missingness
   */
  const VCFValue& justGet(unsigned int i) {
    if (i >= fd.size()) {
      return VCFIndividual::defaultVCFValue;
    }
    return (this->fd[i]);
  };
  
  VCFValue& getSelf() {return this->self;};
  const VCFValue& getSelf() const{
    return this->self;
  };

  size_t size() const {return this->fd.size();}
  /**
   * dump the content of VCFIndividual column
   */
  void output(FILE* fp) const{
    for (unsigned int i = 0; i < fd.size(); ++i){
      if (i)
        fputc(':', fp);
      this->fd[i].output(fp);
    }
  };
private:
  bool inUse;
  std::string name;         // id name
  VCFValue self;            // whole field for the individual (unparsed)
  VCFBuffer parsed;       // store parsed string (where \0 added)
  std::vector<VCFValue> fd; // each field separated by ':'
  static VCFValue defaultVCFValue;
}; // end VCFIndividual

#endif /* _VCFINDIVIDUAL_H_ */
