#ifndef _VCFINFO_H_
#define _VCFINFO_H_

#include <string>
#include "VCFValue.h"
#include "base/OrderedMap.h"

class VCFInfo {
 public:
  VCFInfo() { this->hasParsed = false; };

  void reset() { this->hasParsed = false; };

  /**
   * Check whethere VCF INFO column has @param tag
   * Set @param isMissing to whether this @param tag is missing
   * @return VCFValue of the value of given tag
   * e.g. get("DP", &isMissing) for INFO="DP=3", will return VCFValue of "3" and
   * isMissing = false
   */
  const VCFValue& getTag(const char* tag, bool* isMissing);

  void at(unsigned int idx, std::string* key, VCFValue* value) const {
    if (idx >= this->data.size()) {
    }
    this->data.at(idx, key, value);
  }

  void parse(const VCFValue& v) {
    // this->self = v;
    this->parsed.attach(v.line + v.beg, v.end - v.beg);
    this->hasParsed = false;
  };
  void parseActual() {
    if (this->hasParsed) return;
    // reset data
    this->data.clear();

    // parse key and values
    int state =
        0;  // 0: key, 1: value (indicating current status for value.line[end])
    int end = 0;
    std::string key;
    VCFValue value;
    value.beg = 0;
    value.line = this->parsed.getBuffer();
    const int len = (int)this->parsed.size();
    while (end <= len) {
      if (this->parsed[end] == '=') {
        if (state == 0) {
          key.assign(&(value.line[value.beg]), &(value.line[end]));
          this->parsed[end] = '\0';
          value.beg = end + 1;
          state = 1;
        } else if (state == 1) {
          fprintf(stderr, "Possible wrong format in %s\n",
                  this->parsed.getBuffer());
          return;
        } else {
          fprintf(stderr, "Corrupted state!\n");
          assert(false);
        }
      } else if (this->parsed[end] == ';' || end == len) {
        if (state == 0) {  // key without value: e.g. ;HM3;
          key.assign(&(value.line[value.beg]), &(value.line[end]));
          value.beg = end;
          value.end = end;
          if (key != ".")  // only store non-missing key
            this->data[key] = value;
          value.beg = end + 1;
        } else if (state == 1) {  // key with value: e.g. ;AC=2;
          value.end = end;
          if (key != ".")  // only store non-missing key
            this->data[key] = value;
          value.beg = end + 1;
          state = 0;
        } else {
          fprintf(stderr, "Corrupted state!\n");
          assert(false);
        };

        if (end == len)
          break;
        else
          this->parsed[end] = '\0';
      }

      ++end;
    }
    // finish up
    this->hasParsed = true;
  };
  // const char* toStr() const { return this->self.toStr(); }
  size_t size() {
    this->parseActual();
    return this->data.size();
  }
  void output(FILE* fp, char c) {
    if (this->hasParsed) {
      int n = data.size();
      char sep = ':';
      std::string key;
      VCFValue value;
      for (int i = 0; i < n; ++i) {
        if (i == n - 1) {
          sep = c;
        }
        data.at(i, &key, &value);
        fputs(key.c_str(), fp);
        // check '='
        if (value.beg != value.end) {
          fputc('=', fp);
          value.output(fp, sep);
        } else {
          fputc(sep, fp);
        }
      }
    } else {
      this->parsed.output(fp, c);
    }
  }

 private:
  bool hasParsed;
  VCFBuffer parsed;
  // VCFValue self;
  // std::string parsed;  /// store parsed (where \0 added) string
  OrderedMap<std::string, VCFValue> data;
  const static VCFValue defaultValue;  // Default empty VCFValue
};                                     // VCFInfo

#endif /* _VCFINFO_H_ */
