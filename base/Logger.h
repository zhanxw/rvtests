#ifndef _LOGGER_H_
#define _LOGGER_H_

#include <stdio.h>
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
        if (Logger::fp == NULL) 
            defaultInit();
    };
    Logger(const char* fn, LogLevel l) {
        init(fn, l);
    };
    ~Logger() {
        if (Logger::fp) {
            fclose(Logger::fp);
            Logger::fp = NULL;
        }
        fprintf(stderr, "logger closed\n");
    };

    static int info(const char* fmt, ...);
    static int warn(const char* fmt, ...);
    static int error(const char* fmt, ...);
    static int fatal(const char* fmt, ...);
    static FILE*& getStream();
  public:
    const static LogLevel INFO = _INFO; 
    const static LogLevel WARN = _WARN; 
    const static LogLevel ERROR = _ERROR;
    const static LogLevel FATAL = _FATAL;
  private:
    void defaultInit(){
        init("default.log", Logger::INFO);
    };
    void init(const char* fn, LogLevel l) {
        if (Logger::fp == NULL) {
            Logger::fp = fopen(fn, "w");
            Logger::level = l;
        } else {
            // already have logger.
            if (l != level) {
                info("%s", "Log level changed\n");
                Logger::level = l;
            }
        };
    };
        
  private:
    static FILE* fp;
    static LogLevel level;
};

#endif /* _LOGGER_H_ */
