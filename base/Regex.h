#ifndef _REGEX_H_
#define _REGEX_H_

// We use PCRE here, use 'man pcreposix' for more information
// accordig to http://lh3lh3.users.sourceforge.net/reb.shtml
// PCRE-posix is fast
#include <stdio.h>
#include <pcreposix.h>
#include <string>
#define ERROR_BUF_LEN 64
class Regex {
public:
    /**
     * read pattern like "=Synonymous,=Indel"
     */
    void readPattern(const char* argRegex) {
        int cflags = 0;
        int ret = regcomp(& this->pattern, argRegex, 0);
        if (ret) {
            regerror(ret, & this->pattern, error_buf, ERROR_BUF_LEN);
            fputs(error_buf, stderr);
            exit(1);
        }
        this->initialized = true;
        
    };
    void readPattern(std::string& argRegex) {
        readPattern(argRegex.c_str());
    };
    /**
     */
    bool match(const char* text) {
        if (!this->initialized) return true;
        int eflags = REG_STARTEND;
        int ret = regexec(&this->pattern, text, 0, 0, eflags);
        if (ret == 0) {
            //printf("Match: %s\n", text);
            return true;
        } else if (ret == REG_NOMATCH){
            //printf("Nomatch: %s\n", text);
            return false;
        } else {
            regerror(ret, & this->pattern, error_buf, ERROR_BUF_LEN);
            fputs(error_buf, stderr);
            exit(1);
        }
        return false;
    };

    /**
     * @return if any pattern matches the text, will return true
     */
    bool match(const char* text, int begin, int end) {
        if (!this->initialized) return true;

        size_t nmatch = 1;
        regmatch_t pmatch[nmatch];
        pmatch[0].rm_so = begin;
        pmatch[0].rm_eo = end;
        int eflags = REG_STARTEND;
        int ret = regexec(&this->pattern, text, 0, pmatch, eflags);
        if (ret == 0) {
            //printf("Match: %s\n", text);
            return true;
        } else if (ret == REG_NOMATCH){
            //printf("Nomatch: %s\n", text);
            return false;
        } else {
            regerror(ret, & this->pattern, error_buf, ERROR_BUF_LEN);
            fputs(error_buf, stderr);
            exit(1);
        }
        return false;
    };
    Regex() {
        this->initialized = false;
    }
    ~Regex(){
        if (this->initialized)
            regfree(&pattern);
        this->initialized = false;
    }
private:
    bool initialized;
    regex_t pattern;
    unsigned int numPattern;
    char error_buf[ERROR_BUF_LEN];

};

#endif /* _REGEX_H_ */
