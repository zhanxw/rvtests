#ifndef _VCFINFO_H_
#define _VCFINFO_H_

#include "VCFValue.h"
#include "OrderedMap.h"
#include <string>

class VCFInfo{
public:
  VCFInfo() {
    this->hasParsed = false;
  };

  void reset() { this-> hasParsed = false;};

  /**
   * Check whethere VCF INFO column has @param tag
   * Set @param isMissing to whether this @param tag is missing
   * @return VCFValue of the value of given tag
   * e.g. get("DP", &isMissing) for INFO="DP=3", will return VCFValue of "3" and isMissing = false
   */
  const VCFValue& getTag(const char* tag, bool* isMissing) ;

  void at(unsigned int idx, std::string* key, VCFValue* value) const {
    if (idx >= this->data.size()) {
      
    };
    this->data.at(idx, key, value);
  };
  
  void parse(const VCFValue& v) {
    this->self = v;
    this->hasParsed = false;
  };
  void parseActual(){
    // reset data
    this->data.clear();
    
    // dupliate string
    this->parsed = this->self.line + this->self.beg;
    
    // parse key and values
    int state = 0; // 0: key, 1: value
    unsigned int end = 0;
    std::string key;
    VCFValue value;
    value.beg = 0;
    value.line = this->parsed.c_str();
    while ( end <= this->parsed.size()) {
      if (this->parsed[end] == '=') {
        if (state  == 0) {
          key = this->parsed.substr(value.beg, end - value.beg);
          this->parsed[end] = '\0';
          value.beg = end + 1;
          state = 1;
        } else if (state == 1) {
          fprintf(stderr, "Possible wrong format in %s\n", this->parsed.c_str());
          return;
        } else {
          fprintf(stderr, "Corrupted state!\n");
          assert(false);
        }
      } else if (this->parsed[end] == ';' ||
                 end == this->parsed.size()) {
        if (state == 0) { // key without value: e.g. ;HM3;
          key = this->parsed.substr(value.beg, end - value.beg);
          value.beg = end;
          value.end = end;
          this->data[key] = value;
          value.beg = end + 1;
        } else if (state == 1) { // key with value: e.g. ;AC=2;
            value.end = end;
            this->data[key] = value;
            value.beg = end + 1;
            state = 0;
        } else {
          fprintf(stderr, "Corrupted state!\n");
          assert(false);
        };

        if (end == this->parsed.size())
          break;
        else
          parsed[end] = '\0';
      }

      ++ end;
    }
    // finish up
    this->hasParsed = true;
  };
  const char* toStr() {
    return this->self.toStr();
  }
  size_t size() {
    this->parseActual();
    return this->data.size();
  }
  
private:
  VCFValue self;
  bool hasParsed;
  std::string parsed; /// store parsed (where \0 added) string
  OrderedMap< std::string, VCFValue > data;

  const static VCFValue defaultValue;  // Default empty VCFValue
}; // VCFInfo

#endif /* _VCFINFO_H_ */
