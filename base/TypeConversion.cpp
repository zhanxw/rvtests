#include "TypeConversion.h"
#include "Utils.h"

#include <algorithm>

int chrom2int(const std::string& chrom) {
    int b = 0;
    if (hasLeadingChr(chrom))
        b = 3;
    size_t e;
    e = chrom.find('_', b);
    std::string t = chrom.substr(b, e - b);
    if (t.size() == 0) return -1;
    int ret;
    if (str2int(t.c_str(), &ret)){
        if (e == chrom.npos ){
            return ret;
        } else {
            return (ret + 100);
        }
    } else {
        if ( t == "X" ) return 23;
        if ( t== "Y" ) return 24;
        if ( t== "MT" ) return 25;
        return 1000 + int(t[0]);
    }
}

// convert int to comma-separated string type
// e.g. -123456 => "-123,456"
std::string toStringWithComma(int in){
  std::string ret;
  int plus = in < 0 ? -in : in;
  if (in == INT_MIN) {
    plus = - (1 + in); // note the INT_MIN = -INTMAX -1
  }
  int digits = 0;
  div_t d;
  do {
    d = div(plus, 10);
    plus = d.quot;
    ret.push_back('0' + d.rem);
    digits ++;
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
    ret[0] ++;
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}

