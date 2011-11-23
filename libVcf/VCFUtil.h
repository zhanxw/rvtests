#ifndef _VCFUTIL_H_
#define _VCFUTIL_H_

#include "Exception.h"
#include "RangeList.h"
#include "Utils.h"

#include "VCFConstant.h"
#include "VCFFunction.h"
#include "VCFHeader.h"
#include "VCFInputFile.h"
#include "VCFOutputFile.h"
#include "PlinkOutputFile.h"


#if 0
typedef int(*PARSE_FUNCTION)(const char* s, int beg, int end);
int parseGenotype(const char* line, int beg, int end){
    if (line[beg] == '.') return  MISSING_GENOTYPE;
    int gt = 0;
    if (!isdigit(line[beg])){
        FATAL("Wrong genotype format");
    };
    gt += (line[beg] - '0');
    if (line[beg+1] == ':')  // sex chromosome
        return gt; 
    if (line[beg+1] == '/' || line[beg+1] == '|'){
        if (line[beg+2] == '.') {
            return MISSING_GENOTYPE;
        }
        if (!isdigit(line[beg+2])){
            FATAL("Wrong genotype format");
        };
        gt += (line[beg+2 ] - '0');
        return gt;
    } else {
        FATAL("Wrong genotype format");
    }
    return -1;
};
#endif


#endif /* _VCFUTIL_H_ */
