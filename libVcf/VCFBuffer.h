#ifndef _VCFBUFFER_H_
#define _VCFBUFFER_H_

#include <stdio.h>
#include <stdlib.h>  // for size_t, fprintf
// #include <string.h>  // for strlen, memcpy
#include <cassert>
#include <string>

class FileWriter;

// holding a memory area for a string
// the content may be changed
class VCFBuffer {
 public:
  VCFBuffer() : buf(NULL), len(0), bufLen(0){};
  // VCFBuffer(const char* s) { this->copy(s); };
  void attach(char* s, int l) {
    this->buf = s;
    this->len = this->bufLen = l;
  }
  void attach(std::string& s) { attach(&s[0], (int)s.size()); }
#if 0
  void copy(const char* s) {
    this->len = strlen(s);
    if (this->bufLen < this->len + 1) {
      delete[] this->buf;
      this->buf = NULL;
    }
    if (NULL == this->buf) {
      this->buf = new char[this->len + 1];
      assert(this->buf);
      this->bufLen = this->len + 1;
    }
    memcpy(this->buf, s, this->len);
    this->buf[this->len] = '\0';
  };
  ~VCFBuffer() {
    if (this->buf) {
      delete[] buf;
      this->buf = NULL;
    }
  }
#endif
  char& operator[](const int idx) { return buf[idx]; }
  const char& operator[](const int idx) const { return buf[idx]; }
#if 0
  VCFBuffer& operator=(const char* s) {
    this->copy(s);
    return (*this);
  };
#endif
  void clear() { this->len = 0; };

  const char* c_str() const { return this->buf; };
  char* getBuffer() const { return this->buf; };
  size_t size() const { return this->len; };
  void dump(int firstNumChar) const {
    for (int i = 0; i < firstNumChar; i++) {
      fprintf(stderr, "%d: %c (%d)\n", i, buf[i], buf[i]);
    }
  }
  void output(FILE* fp) const {
    for (size_t i = 0; i != len; ++i) {
      fputc(buf[i], fp);
    }
  }
  void output(FILE* fp, char c) const {
    output(fp);
    fputc(c, fp);
  }
  void output(FileWriter* fp, char c) const;

 private:
  VCFBuffer(VCFBuffer& b);
  VCFBuffer& operator=(const VCFBuffer& b);

 private:
  char* buf;      // pointer to string
  size_t len;     // len (not including the trailing '\0')
  size_t bufLen;  // memory capacity
};

#endif /* _VCFBUFFER_H_ */
