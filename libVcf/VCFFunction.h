#ifndef _VCFFUNCTION_H_
#define _VCFFUNCTION_H_

#include "VCFValue.h"


/* inline int parseTillChar(const char c, const char* line, const int beg, VCFValue* ret) { */
/*     assert(ret); */
/*     ret->line = line; */
/*     ret->beg = beg; */
/*     ret->end = ret->beg; */
/*     while(line[ret->end] != c && line[ret->end] != '\0') { */
/*         ret->end++; */
/*     } */
/*     return ret->end; */
/* }; */

/* inline int parseTillChar(const char* c, const char* line, const int beg, VCFValue* ret) { */
/*     assert(ret); */
/*     ret->line = line; */
/*     ret->beg = beg; */
/*     ret->end = ret->beg; */
/*     int cLen = strlen(c); */
/*     while(!index(c, line[ret->end])){ // line[beg] is not a separator */
/*         ret->end++; */
/*     } */
/*     return ret->end; */
/* }; */

/* /\** */
/*  * @return 0: found @param id in @paraminfo, and store results in @param result */
/*  * @return -1: not found, @param result value is undetermined */
/*  *\/ */
/* inline int parseVCFINFO(std::string* result, const std::string info, const char* id) { */
/*     size_t begin = info.find(id); */
/*     if (begin == std::string::npos) { */
/*         return -1; */
/*     } */
/*     result->clear(); */
/*     assert(info[++begin] == '='); */
/*     size_t len = info.size(); */
/*     while (begin < len && info[begin] != ';') */
/*         (*result).push_back(info[begin++]); */
/* } */

/* /\** */
/*  * @return -1: we cannot find depth */
/*  *\/ */
/* inline int obtainDepth(const char* s, int start, int end) { */
/*     const char* DEPTH_STRING = "DP"; */
/*     const int dpStrLen = strlen(DEPTH_STRING); */
/*     bool isMatch = true; */
/*     int dpBegin, dpEnd; */
/*     for (int i = start; i < end; i++) { */
/*         if (s[i] != DEPTH_STRING[0]) */
/*             continue; */
/*         for (int j = 1; j < dpStrLen; j++) { */
/*             if (s[i+j] != DEPTH_STRING[j]) { */
/*                 isMatch = false; */
/*                 break; */
/*             } */
/*         } */
/*         if (isMatch == false)  */
/*             continue; */
/*         else { */
/*             isMatch = true; */
/*             dpBegin = i; */
/*             break; */
/*         } */
/*     } */
/*     if (!isMatch)  */
/*         return -1; */
    
/*     // move to the first digit */
/*     dpBegin += dpStrLen;  */
/*     assert(s[dpBegin] == '='); */
/*     dpBegin++; */

/*     // put number in buffer and then convert */
/*     static char buffer[32]; // handle 32 digits integer at most */
/*     int i = 0; */
/*     while(s[dpBegin] >= '0' && s[dpBegin] <= '9' && dpBegin < end) { */
/*         buffer[i++] = s[dpBegin++]; */
/*     } */
/*     buffer[i] = '\0'; */
/*     return atoi(buffer); */
/* }; */


#endif /* _VCFFUNCTION_H_ */
