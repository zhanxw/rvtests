#include "TypeConversion.h"
#include "Utils.h"

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
