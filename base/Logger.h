#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <stdio.h>
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

  static int info(const char* fmt, ...);
  static int infoToFile(const char* fmt, ...); // for some logs not shown on the screen, e.g. parameters
  static int warn(const char* fmt, ...);
  static int error(const char* fmt, ...);
  static int fatal(const char* fmt, ...);
  static FILE*& getHandle();
public:
  const static LogLevel INFO = _INFO;
  const static LogLevel WARN = _WARN;
  const static LogLevel ERROR = _ERROR;
  const static LogLevel FATAL = _FATAL;
private:
  void init(const char* fn) {
    if (Logger::fileHandle == NULL) {
      Logger::fileHandle = fopen(fn, "w");
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
