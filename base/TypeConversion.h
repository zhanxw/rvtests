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

// convert double/float to string type
// we try to mimic the '%g' in printf
template<class T>
inline std::string floatToString(T i){
    std::stringstream ss;
    ss.precision(6);
    ss << std::noshowpoint << i;
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

// convert std::string to integer
// @return true if conversion succeed
inline bool str2int(const std::string& input, int* output) {
  return str2int(input.c_str(), output);
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
inline bool str2double(std::string& input, double* output) {
  return str2double(input.c_str(), output);
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
extern int chrom2int(const std::string& chrom);
#endif /* _TYPECONVERSION_H_ */
