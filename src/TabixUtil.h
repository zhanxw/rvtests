#ifndef _TABIXUTIL_H_
#define _TABIXUTIL_H_

#include <string>

int tabixIndexFile(const std::string& fn,
                   int skip = 0,
                   char meta = '#',
                   int chrom = 1,
                   int startPos = 2,
                   int endPos = 0
                   );

#endif /* _TABIXUTIL_H_ */
