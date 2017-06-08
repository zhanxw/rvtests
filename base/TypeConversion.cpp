#include "TypeConversion.h"
#include "Utils.h"

#include <algorithm>

class atoi_func {
 public:
  atoi_func() : value_() {}

  inline int value() const { return value_; }

  inline bool operator()(const char* str, size_t len) {
    value_ = 0;
    int sign = 1;
    if (str[0] == '-') {  // handle negative
      sign = -1;
      ++str;
      --len;
    }

    switch (len) {  // handle up to 10 digits, assume we're 32-bit
      case 10:
        value_ += (str[len - 10] - '0') * 1000000000;
      case 9:
        value_ += (str[len - 9] - '0') * 100000000;
      case 8:
        value_ += (str[len - 8] - '0') * 10000000;
      case 7:
        value_ += (str[len - 7] - '0') * 1000000;
      case 6:
        value_ += (str[len - 6] - '0') * 100000;
      case 5:
        value_ += (str[len - 5] - '0') * 10000;
      case 4:
        value_ += (str[len - 4] - '0') * 1000;
      case 3:
        value_ += (str[len - 3] - '0') * 100;
      case 2:
        value_ += (str[len - 2] - '0') * 10;
      case 1:
        value_ += (str[len - 1] - '0');
        value_ *= sign;
        return value_ > 0;
      default:
        return false;
    }
  }

 private:
  int value_;
};

// convert std::string to integer
// @return true if conversion succeed
bool str2int(const char* input, int* output) {
  // Limitation: this method only works for 32-bit integer
  // We did not use strtol, as it is slower althought it can convert to
  // `long` instead of `int`.
  // The implementaiton may be more compatible than strtol() on musl

  size_t len = 0;
  unsigned char neg = 0;
  unsigned long value = 0;

  if (!input) {
    assert(input && "null input string for str2int()");
    goto err;
  }
  // skip spaces
  while (*input == ' ') {
    input++;
  }

  // check sign
  if (*input == '-') {
    neg = 1;
    input++;
  }

  if (*input == '\0') {
    goto err;
  }

  while (*input == ' ') {
    input++;
  }

  while (input[len] >= '0' && input[len] <= '9') {
    len++;
  }

  // unroll the computation to speed up
  switch (len) {  // handle up to 10 digits, assume we're 32-bit
    case 10:
      value += (input[len - 10] - '0') * 1000000000;
    case 9:
      value += (input[len - 9] - '0') * 100000000;
    case 8:
      value += (input[len - 8] - '0') * 10000000;
    case 7:
      value += (input[len - 7] - '0') * 1000000;
    case 6:
      value += (input[len - 6] - '0') * 100000;
    case 5:
      value += (input[len - 5] - '0') * 10000;
    case 4:
      value += (input[len - 4] - '0') * 1000;
    case 3:
      value += (input[len - 3] - '0') * 100;
    case 2:
      value += (input[len - 2] - '0') * 10;
    case 1:
      value += (input[len - 1] - '0');

      // check range
      // valid 32bit int is from [-2147483648, 2147483647]
      if (neg) {
        if (value > (unsigned long)INT_MAX + 1) {
          goto err;
        } else {
          *output = -value;
          return true;
        }
      } else {
        if (value > (unsigned long)INT_MAX) {
          goto err;
        }
        *output = value;
        return true;
      }

    default:
      return false;
  }
err:
  *output = 0;
  return false;
}

int chrom2int(const std::string& chrom) {
  int b = 0;
  if (hasLeadingChr(chrom)) b = 3;
  size_t e;
  e = chrom.find('_', b);
  std::string t = chrom.substr(b, e - b);
  if (t.size() == 0) return -1;
  int ret;
  if (str2int(t.c_str(), &ret)) {
    if (e == chrom.npos) {
      return ret;
    } else {
      return (ret + 100);
    }
  } else {
    if (t == "X") return 23;
    if (t == "Y") return 24;
    if (t == "MT") return 25;
    return 1000 + int(t[0]);
  }
}

// convert int to comma-separated string type
// e.g. -123456 => "-123,456"
std::string toStringWithComma(int in) {
  std::string ret;
  int plus = in < 0 ? -in : in;
  if (in == INT_MIN) {
    plus = -(1 + in);  // note the INT_MIN = -INTMAX -1
  }
  int digits = 0;
  div_t d;
  do {
    d = div(plus, 10);
    plus = d.quot;
    ret.push_back('0' + d.rem);
    digits++;
    if (digits == 3) {
      ret.push_back(',');
      digits = 0;
    }
    if (d.quot == 0) break;
  } while (plus > 0);
  if (ret[ret.size() - 1] == ',') {
    ret.resize(ret.size() - 1);
  }
  if (in < 0) {
    ret.push_back('-');
  }
  if (in == INT_MIN) {
    ret[0]++;
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}
