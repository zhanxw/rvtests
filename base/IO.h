#ifndef _IO_H_
#define _IO_H_

#define UNUSED(x) ((void)(x))
#include <assert.h>  //assert
#include <stdarg.h>  // va_list
#include <stdio.h>   //fopen
#include <stdlib.h>  //malloc
#include <string.h>  //strchr

#include <string>
#include <vector>

#include "SimpleString.h"

// #define IO_DEBUG

typedef enum FileType {
  PLAIN = 0,
  GZIP = 1,
  BZIP2 = 2,
  BGZIP = 3,
  UNKNOWN = 99
} FileType;

/**
 * Sample usage:
 * std::string line;
 * LineReader* f = AbstractFileReader::open("abc");
 * while (f->readLine(line) > 0) {
 *   ...
 * }
 * delete f;
 */
class AbstractFileReader {
 public:
  typedef enum FileType {
    PLAIN = 0,
    GZIP = 1,
    BZIP2 = 2,
    UNKNOWN = 99
  } FileType;
  virtual ~AbstractFileReader() {
  }  // make it virtual so subclass types can close file handle
  static AbstractFileReader* open(const char* fileName);
  static void close(AbstractFileReader** f);
  // virtual functions
  // each specific file type will need to implement the following function
  // virtual int readLine(std::string* line) = 0;
  // virtual int readLineBySep(std::vector<std::string>* fields, const char*
  // sep) = 0;
  virtual int getc() = 0;
  virtual bool isEof() = 0;
  virtual void close() = 0;
  virtual int read(void* buf, int len) = 0;
  // common utility function
  static FileType checkFileType(const char* fileName);

 protected:
  AbstractFileReader() {}  // forbid explicit create AbstractFileReader class.
};

//////////////////////////////////////////////////////////////////////

/**
 * Example code:

   BufferedReader br("Makefile", 200);
   char buf[500] = {};
   int nRead = 0;
   while ( (nRead = br.read(buf, 500)) > 0) {
     for (int i = 0; i < nRead; i++) {
       printf("%c", buf[i]);
     }
   }
*/
class BufferedReader : public AbstractFileReader {
 public:
  BufferedReader(const char* fileName, int bufferCapacity)
      : bufCap(0), bufEnd(0), bufPtr(0), buf(NULL), fp(NULL), s(1024 * 1024) {
#ifdef IO_DEBUG
    fprintf(stderr, "BufferedReader open %s\n", fileName);
#endif
    // initialize buf
    if (bufferCapacity == 0) {
      fprintf(stderr,
              "Buffer size should be greater than 0, now use default buffer "
              "size 1M instead of %d.\n",
              bufferCapacity);
      this->bufCap = 1024 * 1024;
    } else {
      this->bufCap = (int)(bufferCapacity);
    }
    this->buf = new char[this->bufCap];
    if (!this->buf) {
      fprintf(stderr, "Cannot allocate buffer for BufferedReader.\n");
      return;
    }
    this->bufPtr = 0;
    this->bufEnd = 0;
    // initialize fp
    this->fp = AbstractFileReader::open(fileName);
    if (!this->fp) {
      fprintf(stderr, "Canont open file %s\n", fileName);
      this->fp = NULL;
      // need to quit to prevent further actions
      exit(1);
    }
  }
  virtual ~BufferedReader() { this->close(); }
  int getc() {
    if (this->bufPtr == this->bufEnd) {  // buffer all used, need to refresh
      this->bufEnd = this->fp->read(this->buf, this->bufCap);
      this->bufPtr = 0;
    }

    if (this->bufPtr < this->bufEnd)
      return (this->buf[this->bufPtr++]);
    else
      return EOF;
  }
  bool isEof() {
    // fp reaches the end and read buffer reaches the end
    if (this->fp && this->fp->isEof() && this->bufPtr == this->bufEnd) {
      return true;
    }
    return false;
  }
  void close() {
#ifdef IO_DEBUG
    fprintf(stderr, "BufferedReader close\n");
#endif
    // delete fp
    if (this->fp) {
      AbstractFileReader::close(&fp);
    }
    this->fp = NULL;
    // delete buf
    if (this->buf) {
      delete[] this->buf;
      this->buf = NULL;
      this->bufCap = 0;
      this->bufPtr = 0;
      this->bufEnd = 0;
    }
    this->buf = NULL;
  }
  int read(void* buf, int len) {
    // use current buffer to fill in buf
    int idx = 0;
    while (this->bufPtr < this->bufEnd && len > 0) {
      ((char*)buf)[idx++] = this->buf[this->bufPtr++];
      len--;
    }
    if (len == 0) {
      return idx;
    }
    // fill rest of buf
    int nRead = this->fp->read(((char*)buf) + idx, len);
    idx += nRead;
    // refill buffer
    this->bufEnd = this->fp->read(this->buf, this->bufCap);
    this->bufPtr = 0;
    return idx;
  }
  int readLine(std::string* line) {
    assert(this->fp && line);
    // if (this->fp->isEof()) return 0;
    this->s.resize(0);
    char c;
    unsigned nRead = 0;
    while (true) {
      c = this->getc();
      if (c == EOF) {
        (*line) = s.data();
        return nRead;
      } else if (c == '\r') {
        // skip this
        continue;
      } else if (c == '\n') {
        ++nRead;
        (*line) = s.data();
        return nRead;
      } else {  // normal characters
        ++nRead;
        s.append(c);
      }
    }
    assert(false);  // should not reach here
    return 0;
  }
  int readLineBySep(std::vector<std::string>* fields, const char* sep) {
    assert(this->fp && fields && sep);
    // if (this->fp->isEof()) return 0;
    fields->resize(0);
    char c;
    std::string s;
    s.reserve(4096);  // reserve space
    while (true) {
      c = this->getc();
      if (c == EOF) {
        if (s.size()) {
          fields->push_back(s);
        }
        return fields->size();
      } else if (c == '\r') {
        // skip this
        continue;
      } else if (c == '\n') {
        if (s.size()) {
          fields->push_back(s);
        }
        return fields->size();
      } else if (strchr(sep, c) != NULL) {  // separator
        fields->push_back(s);
        s.resize(0);
      } else {  // normal characters
        s.push_back(c);
      }
    }
    assert(false);  // should not reach here
    return 0;
  }

 private:
  int bufCap;  // capacity of the buffer
  int bufEnd;  // bufPtr should not read beyond bufEnd(incluive)
  int bufPtr;  // from which buffer begins to read
  char* buf;   // [0...bufEnd]
  AbstractFileReader* fp;
  SimpleString s;
};  // end BufferedReader

/** Example code:
 // LineReader lr(fn);
 // while(lr.readLine(&line)>0){
 //         fprintf(stdout, "%s\n", line.c_str());
 // }
 */
class LineReader {
 public:
  LineReader(const std::string& fileName) { init(fileName.c_str()); }
  LineReader(const char* fileName) { init(fileName); }

 private:
  LineReader(const LineReader&);
  LineReader& operator=(const LineReader&);

 public:
  void init(const char* fileName) {
#ifdef IO_DEBUG
    fprintf(stderr, "LineReader open %s\n", fileName);
#endif
    this->fp = new BufferedReader(fileName, 1024 * 1024);
    if (!this->fp) {
      fprintf(stderr, "Canont open file %s\n", fileName);
      this->fp = NULL;
    }
  }
  virtual ~LineReader() {
    if (this->fp) {
      fp->close();
      delete this->fp;
      this->fp = NULL;
    }
#ifdef IO_DEBUG
    fprintf(stderr, "LineReader close\n");
#endif
  }
  // return number of characters read.
  // when reading an empty line, will return 1, as we read '\n', however, line
  // will be empty
  // when reading the end, we will return 0
  int readLine(std::string* line) {
    assert(this->fp && line);
    return this->fp->readLine(line);
  }
  // return number of fields read.
  // when reading an empty line, will return 1, meaning 1 field are read,
  // although its content is empty
  // consecutive separators, e.g. \t\t, will yield empty field
  // when reading to the EOF, will return 0.
  int readLineBySep(std::vector<std::string>* fields, const char* sep) {
    assert(this->fp && fields && sep);
    return this->fp->readLineBySep(fields, sep);
  }

 private:
  BufferedReader* fp;
};

/**
 * @return number of empty elements filtered out
 */
extern int removeEmptyField(std::vector<std::string>* fields);

//////////////////////////////////////////////////////////////////////
// FileWriter related classes
class AbstractFileWriter {
 public:
  /// when open is successful, return 0; else: return non-zero
  virtual int open(const char* fn, bool append = false) = 0;
  virtual void close() = 0;
  virtual int write(const char* s) = 0;
  virtual int writeLine(const char* s) = 0;
  virtual ~AbstractFileWriter() = 0;
};

/**
 * design:
 *  a high level file class, underlying using BufferedFileWriter
 * usage:
 * FileWriter* fout = new FileWriter("a.txt", "w");
 * fout->write("abc");
 * fout->writeLn("abc");
 * fout->close();
 * delete fout->write;
 */
class FileWriter {
 public:
  FileWriter(const char* fileName, bool append = false);
  FileWriter(const char* fileName, FileType t);
  void createBuffer() {
    // create buffer for formatted string
    this->bufLen = 1024;
    this->buf = new char[this->bufLen];
    if (!this->buf) {
      fprintf(stderr, "Cannot allocate printf buffer for FileWriter.\n");
    }
  }
  void close() {
    if (this->fp) {
      this->fp->close();
      delete this->fp;
      this->fp = NULL;
    }
    if (this->fpRaw) {
      delete this->fpRaw;
      this->fpRaw = NULL;
    }
    if (this->buf) {
      delete[] this->buf;
      this->buf = NULL;
    }
#ifdef IO_DEBUG
    fprintf(stderr, "FileWriter desc()\n");
#endif
  }
  ~FileWriter() { this->close(); }
  int write(const char c) {
    char s[] = "0";
    s[0] = c;
    return this->fp->write(s);
  }
  int write(const char* s) { return this->fp->write(s); }
  int write(const std::string& s) { return this->fp->write(s.c_str()); }
  int writeLine(const char* s) {
    int ret = this->fp->write(s);
    this->fp->write("\n");
    return (ret + 1);
  }
  // if @param fileName ends with @param suffix, then return true;
  static bool checkSuffix(const char* fileName, const char* suffix) {
    int lf = strlen(fileName);
    int ls = strlen(suffix);
    if (lf < ls) return false;
    for (int i = lf - ls, j = 0; j < ls;) {
      if (fileName[i++] != suffix[j++]) return false;
    }
    return true;
  }

  /**
   * format string to this->buf, then write it out
   */
  // since in C++, a hidden this pointer is passed, we will use 2, 3 instead of
  // 1 and 2
  int printf(const char* fmt, ...) __attribute__((format(printf, 2, 3))) {
    // we'll put the formatted string  to internal buffer
    va_list args;
    int ret;
    int newBufLen;

    while (1) {
      /* Try to print in the allocated space. */
      va_start(args, fmt);
      ret = vsnprintf(this->buf, this->bufLen, fmt, args);
      va_end(args);
      /* If that worked, return the string. */
      if (ret > -1 && ret < this->bufLen) {
        return this->write(this->buf);
      }
      /* Else try again with more space. */
      if (ret > -1)             /* glibc 2.1 */
        newBufLen = ret + 1;    /* precisely what is needed */
      else                      /* glibc 2.0 */
        newBufLen = bufLen * 2; /* twice the old size */
      increaseBufferTo(newBufLen);
    }
  }

 private:
  void increaseBufferTo(int newBufLen) {
    delete[] this->buf;
    this->buf = new char[newBufLen];
    if (!this->buf) {
      fprintf(stderr, "%s:%d Cannot increase printf buffer for FileWriter.\n",
              __FILE__, __LINE__);
      exit(1);
    }
    this->bufLen = newBufLen;
  }
  AbstractFileWriter* fp;
  AbstractFileWriter* fpRaw;
  char* buf;
  int bufLen;
};  // end class FileWriter

bool fileExists(std::string fn);

#endif /* _IO_H_ */
