#include "Regex.h"

int Regex::readPattern(const char* argRegex) {
  if (this->initialized) {
    regfree(&this->pattern);
    this->initialized = false;
  }
  // int cflags = 0;
  int ret = regcomp(&this->pattern, argRegex, 0);
  if (ret) {
    regerror(ret, &this->pattern, error_buf, ERROR_BUF_LEN);
    fputs(error_buf, stderr);
    return -1;
  }
  this->initialized = true;
  return 0;
};
int Regex::readPattern(const std::string& argRegex) {
  return readPattern(argRegex.c_str());
};
/**
 * @return true if matches.
 */
bool Regex::match(const char* text) {
  if (!this->initialized) {
    fprintf(stderr, "Uninitialized regex!\n");
    return false;
  }
  if (text[0] == '\0') return false;
  int ret = regexec(&this->pattern, text, 1, &this->matchResult, 0);
  if (ret == 0) {
    /* printf("Match: %s\n", text); */
    return true;
  } else if (ret == REG_NOMATCH) {
    // printf("Nomatch: %s\n", text);
    return false;
  } else {
    regerror(ret, &this->pattern, error_buf, ERROR_BUF_LEN);
    fputs(error_buf, stderr);
    return false;
  }
  return false;
};

/**
 * @return if any pattern matches the @param text[begin...end], will return
 * true
 * @param begin: inclusive
 * @param end: exclusive
 */
bool Regex::match(const char* text, int begin, int end) {
  if (!this->initialized) {
    fprintf(stderr, "Uninitialized regex!\n");
    return false;
  }
  if (begin == end) {
    // empty string
    return false;
  };
  /* size_t nmatch = 1; */
  /* regmatch_t pmatch[nmatch]; */
  this->matchResult.rm_so = begin;
  this->matchResult.rm_eo = end;
  int eflags = REG_STARTEND;
  int ret =
      regexec(&this->pattern, text + begin, 1, &this->matchResult, eflags);
  if (ret == 0) {
    // check range
    if (this->matchResult.rm_eo <= end) {
      return true;
    } else {
      return false;
    }
  } else if (ret == REG_NOMATCH) {
    // printf("Nomatch: %s\n", text);
    return false;
  } else {
    regerror(ret, &this->pattern, error_buf, ERROR_BUF_LEN);
    fputs(error_buf, stderr);
    return false;
  }
  return false;
};
Regex::Regex() { this->initialized = false; }
Regex::~Regex() {
  if (this->initialized) regfree(&pattern);
  this->initialized = false;
}
