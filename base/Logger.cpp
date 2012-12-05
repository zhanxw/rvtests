#include "Logger.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

FILE* Logger::fileHandle = NULL;
LogLevel Logger::fileLevel = Logger::INFO;
FILE* Logger::consoleHandle = stderr;
LogLevel Logger::consoleLevel = Logger::INFO;

// use macro to save typing!
#define INFO_LEADING(x)     fputs("[INFO]\t", (x));
#define WARN_LEADING(x)     fputs("[WARN]\t", (x));
#define ERROR_LEADING(x)    fputs("[ERROR]\t", (x));
#define FATAL_LEADING(x)    fputs("[FATAL]\t", (x));
#define PRINT_WITH_FMT(x)   va_list ap;                       \
                            va_start(ap, fmt);                \
                            vfprintf((x), fmt, ap);           \
                            va_end(ap);
#define NEWLINE(x)          fputc('\n', (x));

//////////////////////////////////////////////////////////////////////
// log with format
int Logger::info(const char *fmt, ...) 
{
  if (Logger::fileLevel <= Logger::INFO) {
    INFO_LEADING(Logger::fileHandle);
    PRINT_WITH_FMT(Logger::fileHandle) ;
    NEWLINE(Logger::fileHandle) ;
  }
  if (Logger::consoleLevel <= Logger::INFO) {
    INFO_LEADING(Logger::consoleHandle);
    PRINT_WITH_FMT(Logger::consoleHandle) ;
    NEWLINE(Logger::consoleHandle) ;
  }
  return 0;
}

int Logger::infoToFile(const char *fmt, ...)
{
  if (Logger::fileLevel <= Logger::INFO) {
    INFO_LEADING(Logger::fileHandle);
    PRINT_WITH_FMT(Logger::fileHandle) ;
    NEWLINE(Logger::fileHandle) ;
  }
  return 0;
}

int Logger::warn(const char *fmt, ...)
{
  if (Logger::fileLevel <= Logger::WARN) {
    WARN_LEADING(Logger::fileHandle);
    PRINT_WITH_FMT(Logger::fileHandle) ;
    NEWLINE(Logger::fileHandle) ;
  }
  if (Logger::consoleLevel <= Logger::WARN) {
    WARN_LEADING(Logger::consoleHandle);
    PRINT_WITH_FMT(Logger::consoleHandle) ;
    NEWLINE(Logger::consoleHandle) ;
  }
  return 0;
}
int Logger::error(const char *fmt, ...)
{
  if (Logger::fileLevel <= Logger::ERROR) {
    ERROR_LEADING(Logger::fileHandle);
    PRINT_WITH_FMT(Logger::fileHandle) ;
    NEWLINE(Logger::fileHandle) ;
  }
  if (Logger::consoleLevel <= Logger::ERROR) {
    ERROR_LEADING(Logger::consoleHandle);
    PRINT_WITH_FMT(Logger::consoleHandle) ;
    NEWLINE(Logger::consoleHandle) ;
  }
  return 0;
}
int Logger::fatal(const char *fmt, ...)
{
  if (Logger::fileLevel <= Logger::FATAL) {
    FATAL_LEADING(Logger::fileHandle);
    PRINT_WITH_FMT(Logger::fileHandle) ;
    NEWLINE(Logger::fileHandle) ;
  }
  if (Logger::consoleLevel <= Logger::FATAL) {
    FATAL_LEADING(Logger::consoleHandle);
    PRINT_WITH_FMT(Logger::consoleHandle) ;
    NEWLINE(Logger::consoleHandle) ;
  }
  return 0;
}

FILE*& Logger::getHandle() {
  return Logger::fileHandle;
};
