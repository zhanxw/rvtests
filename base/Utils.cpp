#include "Utils.h"
extern bool hasLeadingChr(const std::string& s) {
    if (s.size() > 3 && 
        (s[0] == 'c' || s[0] == 'C') &&
        (s[1] == 'h' || s[1] == 'H') &&
        (s[2] == 'r' || s[2] == 'R')){
        return true;
    }
    return false;
};
