#ifndef _TYPECONVERSION_H_
#define _TYPECONVERSION_H_

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <limits.h>
#include <math.h> // for HUGE_VALH, HUGE_VALL
#include <sstream>

// convert double/int/byte to string type
template<class T>
inline std::string toString(T i){
    std::stringstream ss;
    ss << i;
    return ss.str();
}

// convert std::string to integer
// @return true if conversion succeed
inline bool str2int(const char* input, int* output) {
    char* endptr;
    long val;
    errno = 0;
    val = strtol(input, &endptr, 10);

    if ((errno == ERANGE && (val == LONG_MAX || val == LONG_MIN))
        || (errno != 0 && val == 0)) {
        perror("strtol");
        return false;
    }

    if (endptr == input) {
        // no digits found
        return false;
    }
    *output = val;
    return true;
}

// convert std::string to double
// @return true if conversion succeed
inline bool str2double(const char* input, double* output) {
    char* endptr;
    double val;

    errno = 0;
    val = strtod(input, &endptr);

    if ((errno == ERANGE && (val == HUGE_VALF || val == HUGE_VALL))
        || (errno != 0 && val == 0.)) {
        perror("strtod");
        return false;
    }

    if (endptr == input) {
        // no digits found
        return false;
    }
    *output = val;
    return true;
}

inline int atoi(const std::string& s) {
    int result;
    bool ret = str2int(s.c_str(), & result);
    if (!ret) {
        return 0;
    }
    return result;
};

inline double atof(const std::string& s) {
    double result;
    bool ret = str2double(s.c_str(), & result);
    if (!ret) {
        return 0.0;
    }
    return result;
};

/**
 * convert @param chrom to integer for easier comparisons
 * leading "chr" will not be considered
 * chr6 -> 6
 * chr6_abc -> 6 + 100
 * chrX -> 23
 * chrY -> 24
 * chrMT -> 25
 * chrOther -> 1000 + ASCII('O')
 * chr -> -1
 */
<<<<<<< HEAD

extern int chrom2int(const std::string& chrom);
=======
extern bool hasLeadingChr(const std::string& );
inline int chrom2int(const std::string& chrom) {
    int b = 0;
    if (hasLeadingChr(chrom))
        b = 3;
    int e;
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
>>>>>>> 4f3ace15dde1c583220b7f36cd802b07a42b3c1c

#endif /* _TYPECONVERSION_H_ */
