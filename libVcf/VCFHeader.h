#ifndef _VCFHEADER_H_
#define _VCFHEADER_H_

#include <vector>
#include "base/Utils.h"

class VCFHeader {
 public:
  /**
   * Use @param s as header
   */
  void setHeader(const std::string& s) {
    stringTokenize(s, "\n", &this->data);
    // filter out emptye entries
    size_t b = 0;
    for (size_t i = 0; i < this->data.size(); ++i) {
      if (this->data[i].empty()) {
        continue;
      }
      this->data[b] = this->data[i];
      b++;
    }
    this->data.resize(b);
    // //debug
    // for (size_t i = 0; i < this->data.size(); ++i) {
    //   fprintf(stderr, "header[%d] = %s\n", (int)i, data[i].c_str());
    // }
  }
  void push_back(const std::string& s) { this->data.push_back(s); }
  void getPeopleName(std::vector<std::string>* p) const {
    if (!p) return;
    if (this->data.size() < 1) return;
    const std::string ln = this->data[this->data.size() - 1];

    std::vector<std::string> fd;
    stringTokenize(ln, "\t", &fd);
    if (fd.size() < 10) return;
    p->clear();
    for (unsigned int i = 9; i < fd.size(); i++) {
      p->push_back(fd[i]);
    }
  };
  int size() const { return this->data.size(); }
  std::string& operator[](int n) { return this->data[n]; };
  const std::string operator[](int n) const { return this->data[n]; };
  std::string at(int n) { return this->data.at(n); };
  const std::string at(int n) const { return this->data.at(n); };
  void clear() { this->data.clear(); };
  void output(FILE* fp) const {
    for (unsigned int i = 0; i < data.size(); ++i) {
      fprintf(fp, "%s\n", data[i].c_str());
    }
  };
  int getPeopleNumber() const {
    const std::string ln = this->data[this->data.size() - 1];
    std::vector<std::string> fd;
    return stringTokenize(ln, "\t", &fd) - 9;
  };

 private:
  std::vector<std::string> data;
};

#endif /* _VCFHEADER_H_ */
