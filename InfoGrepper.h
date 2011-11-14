#ifndef _INFOGREPPER_H_
#define _INFOGREPPER_H_

// We use PCRE here, use 'man pcreposix' for more information
// accordig to http://lh3lh3.users.sourceforge.net/reb.shtml
// PCRE-posix is fast
#include <pcreposix.h>
#define ERROR_BUF_LEN 64
class InfoGrepper {
public:
    /**
     * read pattern like "=Synonymous,=Indel"
     */

    void readPattern(std::string& argInfoGrep) {
        int cflags = 0;
        int ret = regcomp(& this->pattern, argInfoGrep.c_str(), 0);
        if (ret) {
            regerror(ret, & this->pattern, error_buf, ERROR_BUF_LEN);
            fputs(error_buf, stderr);
            exit(1);
        }
        this->initialized = true;
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
    InfoGrepper() {
        this->initialized = false;
    }
    ~InfoGrepper(){
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

#endif /* _INFOGREPPER_H_ */
