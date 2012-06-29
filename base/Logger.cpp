#include "Logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

FILE* Logger::fp = NULL;
LogLevel Logger::level = Logger::INFO;

// use macro to save typing!
#define INFO_LEADING  fputs("[INFO]\t", Logger::fp);
#define WARN_LEADING  fputs("[WARN]\t", Logger::fp);
#define ERROR_LEADING fputs("[ERROR]\t", Logger::fp);
#define FATAL_LEADING fputs("[FATAL]\t", Logger::fp);
#define PRINT_WITH_FMT     va_list ap; \
                           va_start(ap, fmt);                     \
                           int n = vfprintf(Logger::fp, fmt, ap); \
                           va_end(ap);
#define NEWLINE fputc('\n', Logger::fp);

//////////////////////////////////////////////////////////////////////
// log with format
int Logger::info(const char *fmt, ...)
{
    if (Logger::level > Logger::INFO) return 0;
    INFO_LEADING ;
    PRINT_WITH_FMT ;
    NEWLINE ;
    return n;
}

int Logger::warn(const char *fmt, ...)
{
    if (Logger::level > Logger::WARN) return 0;
    WARN_LEADING ;
    PRINT_WITH_FMT ;
    NEWLINE ;
    return n;
}
int Logger::error(const char *fmt, ...)
{
    if (Logger::level > Logger::ERROR) return 0;
    ERROR_LEADING ;
    PRINT_WITH_FMT ;
    NEWLINE ;
    return n;
}
int Logger::fatal(const char *fmt, ...)
{
    if (Logger::level > Logger::FATAL) return 0;
    FATAL_LEADING ;
    PRINT_WITH_FMT ;
    NEWLINE ;
    return n;
}

FILE*& Logger::getStream() {
    return Logger::fp;
};
