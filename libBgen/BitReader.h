#ifndef _BITREADER_H_
#define _BITREADER_H_

class BitReader {
 public:
  BitReader(uint8_t* ptr, int len, const int B)
      : ptr(ptr), offset(0), B(B), len(len) {
    assert(1 <= B && B <= 32);
    scale = 1.0;
    for (int i = 0; i < B; ++i) {
      scale *= 2;
    }
    scale -= 1;
    scale = 1.0 / scale;
  }
  // idx is 0 based
  uint32_t getBit(int idx) {
    assert(0 <= idx && idx < len * 8);
    // printf("idx = %d offset = %d len = %d\n", idx, offset, len);
    int byte = (offset + idx) / 8;
    int remain = (offset + idx) % 8;
    return (ptr[byte] & (1 << remain)) != 0 ? 1 : 0;
  }
  float next() {
    int res = 0;
    for (int i = 0; i < B; ++i) {
      res |= getBit(i) << i;
    }
    offset += B;
    return ((float)res) * scale;
  }

 private:
  uint8_t* ptr;
  int offset;
  const int B;
  float scale;
  int len;
};

#endif /* _BITREADER_H_ */
