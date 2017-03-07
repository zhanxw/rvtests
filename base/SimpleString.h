#ifndef _SIMPLESTRING_H_
#define _SIMPLESTRING_H_

#include <string.h>  // memcpy

#ifndef NDEBUG
#include <stdio.h>
#include <stdlib.h>
#endif

class SimpleString {
 public:
  SimpleString(int l) {
    beg = new char[l + 1];
    beg[l] = '\0';
    current = beg;
    cap = l;
  }
  ~SimpleString() { delete[] beg; }
  int reserve(int l) {
    if (l <= 0) return 0;
#ifndef NDEBUG
    if (l > cap) {
      fprintf(stderr, "%s:%d reallocate memory\n", __FILE__, __LINE__);
    }
#endif
    int newCap = l * 2;
    char* newBeg = new char[newCap];
    newBeg[newCap - 1] = '\0';
    memcpy(newBeg, beg, current - beg);
    cap = newCap;
    current = newBeg + (current - beg);
    delete[] beg;
    beg = newBeg;

    return 0;
  }
  int reserveDouble() { return reserve(2 * cap); }
  int resize(int l) {
    if (l > cap) {
#ifndef NDEBUG
      fprintf(stderr, "need reallocation!\n");
#endif
      reserve(l);
    }
    current = beg + l;
    return 0;
  }
  void append(char c) {
    if (current >= beg + cap) {
      reserveDouble();
    }
    *(current) = c;
    current++;
  }
  // append buf[start...(end-1)] to the buf
  void append(char* buf, int start, int end) {
    if (start <= end) {
      return;
    }
    while (current + (end - start) >= beg + cap) {
      reserveDouble();
    }
    const int len = (end - start);
    memcpy(current, buf + start, len);
    current += len;
  }
  const char* data() {
    *current = '\0';
    return this->beg;
  }
  size_t size() const { return current - beg; }
  bool empty() const { return current == beg; }

 private:
  unsigned int next_pow2(unsigned int x) {
    x -= 1;
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);

    return x + 1;
  }

  SimpleString(const SimpleString&);
  SimpleString& operator=(const SimpleString&);

 private:
  char* beg;      // pointer to the data
  char* current;  // current read/write position
  // char* end;      // boundary of the read/write limit
  int cap;  // memory capacity
};

#endif /* _SIMPLESTRING_H_ */
