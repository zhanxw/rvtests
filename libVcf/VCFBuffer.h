#ifndef _VCFBUFFER_H_
#define _VCFBUFFER_H_

#include <cassert>

class VCFBuffer {
public:
VCFBuffer():buf(NULL), len(0), bufLen(0) {};
  VCFBuffer(const char* s) {
    this->copy(s);
  };
  void copy(const char* s){
    this->len = strlen(s);
    if (this->bufLen < this->len + 1) {
      delete[] this->buf;
      this->buf = NULL;
    }
    if (NULL == this->buf){
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
  char& operator[] (const int idx) {
    return buf[idx];
  }
  const char& operator[] (const int idx) const {
    return buf[idx];
  }
  VCFBuffer& operator=(const char* s) {
    this->copy(s);
    return (*this);
  };
  void clear() {
    this->len = 0;
  };
  const char* c_str() const {return  this->buf;};
  size_t size() const {return this->len;};
  void dump(int firstNumChar) const {
    for (int i = 0; i < firstNumChar; i++) {
      fprintf(stderr, "%d: %c (%d)\n", i, buf[i], buf[i]);
    }
  }
private:
  VCFBuffer(VCFBuffer& b);
  VCFBuffer& operator=(const VCFBuffer& b);
  
private:
  char* buf;       // pointer to string
  size_t len;      // len (not including the trailing '\0'
  size_t bufLen;   // memory capacity
};

#endif /* _VCFBUFFER_H_ */
