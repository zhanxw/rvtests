#ifndef _EXCEPTION_H_
#define _EXCEPTION_H_

#include <stdio.h>
#include <stdlib.h>

inline void REPORT(const char* x) { 
    fprintf(stderr, "Report '%s' to zhanxw@umich.edu\n", x ); 
};

inline void FATAL(const char* x) {
    REPORT(x);
    abort();
};

#endif /* _EXCEPTION_H_ */
