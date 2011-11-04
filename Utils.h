#ifndef _UTILS_H_
#define _UTILS_H_

#include <string.h> // for strlen
#include <stdlib.h>

#include <string>
#include <cassert>
#include <algorithm>
/**
 * @return 0: found @param id in @paraminfo, and store results in @param result
 * @return -1: not found, @param result value is undetermined
 */
inline int parseVCFINFO(std::string* result, const std::string info, const char* id) {
    size_t begin = info.find(id);
    if (begin == std::string::npos) {
        return -1;
    }
    result->clear();
    assert(info[++begin] == '=');
    size_t len = info.size();
    while (begin < len && info[begin] != ';')
        (*result).push_back(info[begin++]);
}

/**
 * @return -1: we cannot find depth
 */
int obtainDepth(const char* s, int start, int end) {
    const char* DEPTH_STRING = "DP";
    const int dpStrLen = strlen(DEPTH_STRING);
    bool isMatch = true;
    int dpBegin, dpEnd;
    for (int i = start; i < end; i++) {
        if (s[i] != DEPTH_STRING[0])
            continue;
        for (int j = 1; j < dpStrLen; j++) {
            if (s[i+j] != DEPTH_STRING[j]) {
                isMatch = false;
                break;
            }
        }
        if (isMatch == false) 
            continue;
        else {
            isMatch = true;
            dpBegin = i;
            break;
        }
    }
    if (!isMatch) 
        return -1;
    
    // move to the first digit
    dpBegin += dpStrLen; 
    assert(s[dpBegin] == '=');
    dpBegin++;

    // put number in buffer and then convert
    static char buffer[32]; // handle 32 digits integer at most
    int i = 0;
    while(s[dpBegin] >= '0' && s[dpBegin] <= '9' && dpBegin < end) {
        buffer[i++] = s[dpBegin++];
    }
    buffer[i] = '\0';
    return atoi(buffer);
};

// transpose matrix with dimension @param nr by @param nc
template<class T>
void transposeMatrix(T* m, int nr, int nc) {
    for(int i = 0; i < nr; i++) {
        for (int j = 0; j < i; j++ ){
            std::swap( (*m)[i*nc+j], (*m)[j*nc+i]);
        }
    }
}

// remove the leading 'chr' if any
std::string chopChr(const std::string& s) {
    if (s.size() > 3 && 
        (s[0] == 'c' || s[0] == 'C') &&
        (s[1] == 'h' || s[1] == 'H') &&
        (s[2] == 'r' || s[2] == 'R')){
        return s.substr(3);
    }
    return s;
};

// remove the leading and trailing white spaces 
std::string stringStrip(const std::string& s){
    unsigned int beg = s.find_first_not_of(' ');
    unsigned int end = s.find_last_not_of(' ');
    return s.substr(beg, end-beg);
}

/** tokenize the string
    @return number of tokens we obtained
    e.g. For empty input string, we will return 1, and result will have only 1 element (the empty string)
*/
int stringTokenize(const std::string& str, const std::string& delim, std::vector<std::string>* result){
    assert(result);
    result->clear();
    std::string s;
    unsigned int l = str.size();
    unsigned int i = 0;
    while (i < l) {
        if (delim.find(str[i]) != std::string::npos) { // it's a delimeter
            result->push_back(s);
            s.clear();
        } else {
            s.push_back(str[i]);
        }
        ++i;
    };
    result->push_back(s);    
    return result->size();
};
int stringTokenize(const std::string& str, const char delim, std::vector<std::string>* result){
    std::string d;
    d.push_back(delim);
    return (stringTokenize(str, d, result));
};

/**
 * print out the content for debug only
 */
void dumpStringVector(const std::vector<std::string> s) {
    for (unsigned int i = 0; i < s.size(); i++) {
        fprintf(stdout, "%u: %s\n", i, s[i].c_str());
    }
};
#endif /* _UTILS_H_ */
