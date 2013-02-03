#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

/* class LogLevel{ */
/*   public: */
/*     const static int INFO = 0; */
/*     const static int WARN = 1; */
/*     const static int ERROR = 2; */
/*     const static int FATAL = 3; */
/* }; */
typedef enum{
  _INFO = 0,
  _WARN = 1,
  _ERROR = 2,
  _FATAL = 3
} LogLevel;

class Logger{
public:
  Logger() {
    if (Logger::fileHandle == NULL) {
      char buf [100];
      sprintf(buf, "default.%d.log", (int)(getpid()));
      init(buf);
    }
  };
  Logger(const char* fn) {
    init(fn);
  };
  ~Logger() {
    if (Logger::fileHandle) {
      fclose(Logger::fileHandle);
      Logger::fileHandle = NULL;
    }
  };

  static int info(const char* fmt, ...)   __attribute__ ((format (printf, 1, 2)));
  static int infoToFile(const char* fmt, ...) __attribute__ ((format (printf, 1, 2))); // for some logs not shown on the screen, e.g. parameters
  static int warn(const char* fmt, ...) __attribute__ ((format (printf, 1, 2)));
  static int error(const char* fmt, ...) __attribute__ ((format (printf, 1, 2)));
  static int fatal(const char* fmt, ...) __attribute__ ((format (printf, 1, 2)));
  static FILE*& getHandle();
  static void sync(FILE* f) {
    if (f) {
      fflush(f);
      int no = (fileno(f));
      if ( no > 0) {
        fsync(no);
      }
    }
  }
  static void sync() {
    sync(fileHandle);
    sync(consoleHandle);
  }
public:
  const static LogLevel INFO = _INFO;
  const static LogLevel WARN = _WARN;
  const static LogLevel ERROR = _ERROR;
  const static LogLevel FATAL = _FATAL;
private:
  void init(const char* fn) {
    if (Logger::fileHandle == NULL) {
      Logger::fileHandle = fopen(fn, "w");
      if (Logger::fileHandle == NULL) {
        fprintf(stderr, "Cannot open log file [ %s ].\n", fn);
        exit(1);
      }
    }
    Logger::consoleHandle = stderr;
    Logger::fileLevel = Logger::INFO;
    Logger::consoleLevel = Logger::INFO;
  };
  void setConsoleLevel(LogLevel l){
    if (l != Logger::consoleLevel && l != Logger::consoleLevel) {
      info("%s", "Log console level changed\n");
    }
    Logger::consoleLevel = l;
  };
  void setFileLevel(LogLevel l){
    if (l != Logger::fileLevel && l != Logger::fileLevel) {
      info("%s", "Log file level changed\n");
    }
    Logger::fileLevel = l;
  };
private:
  static FILE* fileHandle;
  static LogLevel fileLevel;  
  static FILE* consoleHandle;
  static LogLevel consoleLevel;
};

#endif /* _LOGGER_H_ */
