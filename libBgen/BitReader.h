#ifndef _BITREADER_H_
#define _BITREADER_H_

/**
 * B = 3
 * data[0]:   22 111 000
 * data[1]:   5 444 333 2
 * data[2]:   777 666 55
 */
class BitReader {
 public:
  BitReader(uint8_t* data, int len, const int B)
      : data(data), offset(0), len(len), availableBits(0), B(B) {
    assert(1 <= B && B <= 32);
    value = 0;
    mask = (1 << B) - 1;
    scale = 1.0;
    for (int i = 0; i < B; ++i) {
      scale *= 2;
    }
    scale -= 1;
    scale = 1.0 / scale;
  }
#if 0
  // idx is 0 based
  uint32_t getBit(int idx) {
    assert(0 <= idx && idx < len * 8);
    // printf("idx = %d offset = %d len = %d\n", idx, offset, len);
    int byte = (offset + idx) / 8;
    int remain = (offset + idx) % 8;
    return (data[byte] & (1 << remain)) != 0 ? 1 : 0;
  }
#endif
  float next() {
    if (B == 8) {
      return (float(data[offset++]) * scale);
    }
    if (B == 16) {
      uint16_t v = *(uint16_t*)(data + offset);
      offset += 2;
      return (float(v) * scale);
    }
    if (B == 32) {
      uint32_t v = *(uint32_t*)(data + offset);
      offset += 4;
      return (float(v) * scale);
    }

    // bits layout:
    // value =         |  0...0   |  to_be_used_bits_1  |
    // data[offset] =  |  ....    |  to_be_used_bits_2  |
    // load data to the left of value:
    // value =    |  ....    |  to_be_used_bits_2  |  to_be_used_bits_1  |
    while (availableBits < B && offset < len) {
      value |= ((uint64_t)data[offset]) << availableBits;
      offset++;
      availableBits += 8;
    }
    float res = float(value & mask);
    availableBits -= B;
    value >>= B;
    return (res * scale);
  }

 private:
  uint8_t* data;
  int offset;
  int len;  // length of the data
  uint8_t availableBits;
  const int B;  // number of bits representing a probability
  uint64_t value;
  uint64_t mask;
  float scale;
};

#endif /* _BITREADER_H_ */
