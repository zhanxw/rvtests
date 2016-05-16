#ifndef _REGEX_H_
#define _REGEX_H_

// We use PCRE here, use 'man pcreposix' for more information
// accordig to http://lh3lh3.users.sourceforge.net/reb.shtml
// PCRE-posix is fast
#include <stdio.h>
#include <string>

// use third/ pcre library instead of
// system provided header file
#include "pcreposix.h"

#define ERROR_BUF_LEN 64
class Regex {
 public:
  /**
   * read pattern like "=Synonymous,=Indel"
   */
  int readPattern(const char* argRegex);
  int readPattern(const std::string& argRegex);
  /**
   * @return true if matches.
   */
  bool match(const char* text);

  /**
   * @return if any pattern matches the @param text[begin...end], will return
   * true
   * @param begin: inclusive
   * @param end: exclusive
   */
  bool match(const char* text, int begin, int end);
  Regex();
  ~Regex();
  bool isInitialized() const { return this->initialized; }

 private:
  bool initialized;
  regex_t pattern;
  char error_buf[ERROR_BUF_LEN];
  regmatch_t matchResult;
};

#endif /* _REGEX_H_ */
