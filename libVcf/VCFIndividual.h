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

    if (fd.size() == 0) {
      fprintf(stderr, "Empty individual column - very strange!!");
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
  VCFValue& operator [] (const unsigned int i) {
    if (i >= fd.size()){
      FATAL("index out of bound!");
    }
    return (this->fd[i]);
  };
  VCFValue& getSelf() {return this->self;};
  void rebuildString(std::string* s) {
    assert(s);
    s->clear();
    for (int i = 0; i < fd.size(); i++){
      if (i)
        s->push_back(':');
      *s += this->fd[i].toStr();
    }
  };
  const VCFValue& getSelf() const{
    return this->self;
  };
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
}; // end VCFIndividual

#endif /* _VCFINDIVIDUAL_H_ */
