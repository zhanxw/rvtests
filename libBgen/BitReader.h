#ifndef _BITREADER_H_
#define _BITREADER_H_

/**
 * B = 3
 * ptr[0]:   22 111 000
 * ptr[1]:   5 444 333 2
 * ptr[2]:   777 666 55
 */
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
    if (B == 8) {
      offset += B;
      return (float(ptr[offset / 8]) * scale);
    }
    if (B == 16) {
      uint16_t v = *(uint16_t*)(ptr + (offset / 8));
      offset += B;
      return (float(v) * scale);
    }
    int res = 0;
    for (int i = 0; i < B; ++i) {
      res |= getBit(i) << i;
    }
    offset += B;
    return ((float)res) * scale;
  }

 private:
  uint8_t* ptr;
  int offset;  // 0-based index
  const int B;
  float scale;
  int len;  // length of the ptr
};

#endif /* _BITREADER_H_ */
