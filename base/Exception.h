#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_

#include <stdio.h>
#include <stdlib.h>

void REPORT(const char* x) { 
    fprintf(stderr, "Report '%s' to zhanxw@umich.edu\n", x ); 
};

void FATAL(const char* x) {
    REPORT(x);
    abort();
};

#endif /* _EXCEPTION_H_ */
