#include "IO.h"

// cannot forward declare an typdef anonymous struct
// http://stackoverflow.com/questions/804894/forward-declaration-of-a-typedef-in-c
// so include the header file
#include "third/samtools/bgzf.h"

#include <algorithm>
#include "base/Utils.h"

//////////////////////////////////////////////////
// Plain file reader
class PlainFileReader : public AbstractFileReader {
 public:
  PlainFileReader(const char* fileName) : fp(NULL) {
    this->open(fileName);
#ifdef IO_DEBUG
    fprintf(stderr, "PlainFileReader() open %s\n", fileName);
#endif
  }
  virtual ~PlainFileReader() {
#ifdef IO_DEBUG
    fprintf(stderr, "~PlainFileReader() close\n");
#endif
    this->close();
  }

  // get a char, if EOF, return EOF
  int getc() { return ::getc(this->fp); }
  // check eof
  bool isEof() { return (feof(this->fp) != 0); }
  // open
  FILE* open(const char* fileName) {
    this->fp = fopen(fileName, "r");
    if (!this->fp) {
      fprintf(stderr, "ERROR: Cannot open file %s\n", fileName);
      this->close();
      exit(1);
    }
    return this->fp;
  }
  // close
  void close() {
    if (this->fp) {
      fclose(fp);
      fp = NULL;
    }
  }
  int read(void* buf, int len) {
    return ::fread(buf, sizeof(char), len, this->fp);
  }

 private:
  FILE* fp;
};

//////////////////////////////////////////////////////////////////////
// Gzip reading class
#include <zlib.h>
class GzipFileReader : public AbstractFileReader {
 public:
  GzipFileReader(const char* fileName) : fp(NULL) {
    this->open(fileName);
#ifdef IO_DEBUG
    fprintf(stderr, "GzipFileReader() open %s\n", fileName);
#endif
  }
  virtual ~GzipFileReader() {
#ifdef IO_DEBUG
    fprintf(stderr, "~GzipFileReader() close\n");
#endif
    this->close();
  }

  // get a char, if EOF, return EOF
  int getc() { return gzgetc(this->fp); }
  // check eof
  bool isEof() { return (gzeof(this->fp) != 0); }
  // open
  gzFile open(const char* fileName) {
    this->fp = gzopen(fileName, "r");
    if (!this->fp) {
      fprintf(stderr, "ERROR: Cannot open gzip file %s\n", fileName);
      this->close();
      exit(1);
    }
    return this->fp;
  }
  // close
  void close() {
    if (this->fp) {
      gzclose(fp);
      fp = NULL;
    }
  }
  int read(void* buf, int len) { return gzread(this->fp, buf, len); }

 private:
  gzFile fp;
};

//////////////////////////////////////////////////////////////////////
// Bzip2 reading class
#include <bzlib.h>
class Bzip2FileReader : public AbstractFileReader {
 public:
  Bzip2FileReader(const char* fileName) : fp(NULL) {
    this->open(fileName);
#ifdef IO_DEBUG
    fprintf(stderr, "Bzip2FileReader() open %s\n", fileName);
#endif
  }
  virtual ~Bzip2FileReader() {
#ifdef IO_DEBUG
    fprintf(stderr, "~Bzip2FileReader() close\n");
#endif
    if (this->fp) {
      BZ2_bzclose(fp);
    }
  }

  // get a char, if EOF, return EOF
  int getc() {
    char c;
    if (this->bzerror != BZ_STREAM_END) {
      int nBuf = BZ2_bzRead(&this->bzerror, this->bzp, &c, sizeof(char));
      if (nBuf) {
        fprintf(stderr, "Read %c \n", c);
        return c;
      }
    }
    return EOF;
  }
  // check eof
  bool isEof() { return (this->bzerror == BZ_STREAM_END); }
  // open
  BZFILE* open(const char* fileName) {
    this->fp = fopen(fileName, "rb");
    if (!this->fp) {
      fprintf(stderr, "ERROR: Cannot open file %s\n", fileName);
      this->close();
      exit(1);  // return NULL;
    }
    this->bzp = BZ2_bzReadOpen(&this->bzerror, this->fp, 0, 0, NULL, 0);

    if (this->bzerror != BZ_OK) {
      BZ2_bzReadClose(&bzerror, this->bzp);
      fprintf(stderr, "ERROR: Cannot open bzip2 file %s\n", fileName);
      this->close();
      return NULL;
    }
    return this->bzp;
  }
  // close
  void close() {
    if (this->bzerror != BZ_STREAM_END) {
      BZ2_bzReadClose(&this->bzerror, this->bzp);
      /* omit code to handle error */
    } else {
      BZ2_bzReadClose(&this->bzerror, this->bzp);
    }
    if (this->fp) fclose(this->fp);
    this->fp = NULL;
    this->bzp = NULL;
    this->bzerror = 0;
  }
  int read(void* buf, int len) {
    return BZ2_bzRead(&this->bzerror, this->bzp, buf, len);
  }

 private:
  FILE* fp;
  BZFILE* bzp;
  int bzerror;
};

#ifdef _USE_KNETFILE
#include "knetfile.h"

class KnetFileReader : public AbstractFileReader {
 public:
  KnetFileReader(const char* fileName) : bgzf_fp(NULL), knet_fp(NULL) {
    this->open(fileName);
#ifdef IO_DEBUG
    fprintf(stderr, "KnetFileReader() open %s\n", fileName);
#endif
  };
  virtual ~KnetFileReader() {
#ifdef IO_DEBUG
    fprintf(stderr, "~KnetFileReader() close\n");
#endif
    this->close();
  };

  // get a char, if EOF, return EOF
  int getc() {
    if (bgzfMode) return bgzf_getc(this->bgzf_fp);
    char c;
    knet_read(this->knet_fp, &c, 1);
    return c;
  }
  // check eof
  bool isEof() {
    if (bgzfMode) {
      if (!this->bgzf_fp) {
        return true;
      }
    } else {
      if (!knet_fp) {
        return true;
      }
    }

    // this is always false, as we don't know the exact file size
    return false;

    // bgzf_check_EOF will 'fseek' to the file end, and check the last 28 bytes
    // but for net resources, it's not possible to reach file end.
    // return bgzf_check_EOF(this->fp);
  }
  // open
  void* open(const char* fileName) {
    size_t l = strlen(fileName);
    if (l > 3 && !strcmp(fileName + l - 3, ".gz")) {
      bgzfMode = true;
    } else {
      bgzfMode = false;
    }
    if (bgzfMode) {
      this->bgzf_fp = bgzf_open(fileName, "r");
      if (!this->bgzf_fp) {
        fprintf(stderr, "ERROR: Cannot open knetfile in bgzf mode: %s\n",
                fileName);
        this->close();
        exit(1);
      }
      return this->bgzf_fp;
    }
    this->knet_fp = knet_open(fileName, "r");
    if (!this->knet_fp) {
      fprintf(stderr, "ERROR: Cannot open knetfile: %s\n", fileName);
      this->close();
      exit(1);
    }
    return this->knet_fp;
  }
  // close
  void close() {
    if (this->bgzfMode) {
      if (this->bgzf_fp) {
        bgzf_close(this->bgzf_fp);
        bgzf_fp = NULL;
      }
    }
    if (this->knet_fp) {
      knet_close(this->knet_fp);
      knet_fp = NULL;
    }
  }
  int read(void* buf, int len) {
    if (bgzfMode) return bgzf_read(this->bgzf_fp, buf, len);
    return knet_read(this->knet_fp, buf, len);
  }

 private:
  bool bgzfMode;
  BGZF* bgzf_fp;
  knetFile* knet_fp;
};  // end KnetFileReader
#endif

//////////////////////////////////////////////////
// BufferedReader
//////////////////////////////////////////////////
BufferedReader::BufferedReader(const char* fileName, int bufferCapacity)
    : bufCap(0), bufEnd(0), bufPtr(0), buf(NULL), fp(NULL) {
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
    fprintf(stderr, "Cannot allocate buffer for BufferedReader. - Exit!\n");
    exit(1);
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

void BufferedReader::close() {
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

bool BufferedReader::isEof() {
  // fp reaches the end and read buffer reaches the end
  if (this->bufPtr < this->bufEnd) return false;

  if (this->fp && this->fp->isEof()) {
    return true;
  }
  return false;
}

int BufferedReader::getc() {
  if (this->bufPtr < this->bufEnd) return (this->buf[this->bufPtr++]);

  // buffer all used, need to refresh
  if (this->bufPtr == this->bufEnd) {
    this->bufEnd = this->fp->read(this->buf, this->bufCap);
    this->bufPtr = 0;
  }
  if (this->bufPtr < this->bufEnd)
    return (this->buf[this->bufPtr++]);
  else
    return EOF;
}

void BufferedReader::refill() {
  bufPtr = bufEnd = 0;
  bufEnd = this->fp->read(this->buf, this->bufCap);
  if (bufEnd < 0) {  // error happens
    bufEnd = 0;
  }
}
/**
 *  0     1     2      3            (4)
 *        ^bufPtr      ^bufEnd      ^bufCap
 *
 * buf[0..3]: data area
 * buf[1] = buf[butPtr]: next position to read data
 * buf[3] = buf[bufEnd]:
 * only buf[0..2] are available
 * if len <=
 * if len > (4 -1), then only raed
 */
int BufferedReader::read(void* pDest, int destLen) {
  if (destLen <= 0) return 0;

  int nRead = 0;
  char* dest = (char*)pDest;

  // process buf[bufPtr..bufEnd)
  int availableData = bufEnd - bufPtr;
  if (availableData >= destLen) {
    memcpy(dest, buf + bufPtr, destLen);
    bufPtr += destLen;
    return destLen;
  }

  // copy all existing data (buf[bufPtr..bufEnd])
  memcpy(dest, buf + bufPtr, availableData);
  nRead += availableData;
  bufPtr += availableData;
  dest += availableData;
  destLen -= availableData;
  assert(bufPtr == bufEnd);

  // refill new data and then copy to @param dest
  while (true) {
    refill();
    if (bufEnd == 0) {  // file end
      return nRead;
    }
    if (bufEnd >= destLen) {
      memcpy(dest, buf, destLen);
      bufPtr += destLen;
      nRead += destLen;
      return nRead;
    }
    memcpy(dest, buf, bufEnd);
    nRead += bufEnd;
    bufPtr = bufEnd;
    dest += bufEnd;
    destLen -= bufEnd;
    assert(bufPtr == bufEnd);
  }

#if 0
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
#endif
}

int BufferedReader::search(int left, int right, const char* sep) {
  assert(right <= bufEnd);
  const char* p;
  for (int i = left; i < right; ++i) {
    p = ssechr(sep, buf[i]);
    if (p != NULL) {
      bufPtr = i + 1;
      return i;
    }
  }
  bufPtr = bufEnd;
  assert(right == bufEnd);
  return bufPtr;
}

int BufferedReader::search(int left, int right, const char* sep1,
                           const char* sep2) {
  assert(right <= bufEnd);
  const char* p;
  for (int i = left; i < right; ++i) {
    p = ssechr(sep1, buf[i]);
    if (p != NULL) {
      bufPtr = i + 1;
      return i;
    }
    p = ssechr(sep2, buf[i]);
    if (p != NULL) {
      bufPtr = i + 1;
      return i;
    }
  }
  bufPtr = bufEnd;
  assert(right == bufEnd);
  return bufPtr;
}

int BufferedReader::readLine(std::string* line) {
  assert(this->fp && line);

  line->reserve(2048);
  line->resize(0);

  int oldPtr, ptr;
  while (true) {
    oldPtr = bufPtr;
    ptr = search(bufPtr, bufEnd, "\r\n");
    line->append(buf + oldPtr, ptr - oldPtr);
    if (ptr == bufEnd) {  // not found
      refill();
      if (bufEnd == 0) {  // file end
        return line->size();
      }
    } else if (buf[ptr] == '\r') {
    } else if (buf[ptr] == '\n') {
      return line->size();
    }
  }
  assert(false);
  return 0;
}

int BufferedReader::readLineBySep(std::vector<std::string>* fields,
                                  const char* sep) {
  assert(this->fp && fields && sep);
  fields->reserve(128);
  fields->resize(1);
  fields->back().resize(0);

  int oldPtr, ptr;
  while (true) {
    oldPtr = bufPtr;
    ptr = search(bufPtr, bufEnd, sep, "\r\n");
    std::string& field = fields->back();

    field.append(buf + oldPtr, ptr - oldPtr);
    if (ptr == bufEnd) {  // not found
      refill();
      if (bufEnd == 0) {               // reach file end
        if (fields->back().empty()) {  // just appended an empty element
          fields->resize(fields->size() - 1);
        }
        return fields->size();
      }
    } else if (ssechr(sep, buf[ptr])) {  // separator
      fields->resize(fields->size() + 1);
      fields->back().resize(0);
    } else if (buf[ptr] == '\r') {
    } else if (buf[ptr] == '\n') {
      return fields->size();
    }
  }
}

//////////////////////////////////////////////////
class TextFileWriter : public AbstractFileWriter {
 public:
  TextFileWriter(const std::string& fn, bool append = false) {
    if (this->open(fn.c_str(), append)) {
      fprintf(stderr, "Cannot create text file %s\n", fn.c_str());
    }
  }
  virtual ~TextFileWriter() {
#ifdef IO_DEBUG
    fprintf(stderr, "TextFileWriter desc()\n");
#endif
    this->close();
  }
  int open(const char* fn, bool append = false) {
    if (append)
      this->fp = fopen(fn, "a");
    else
      this->fp = fopen(fn, "w");
    if (!this->fp) {
      fprintf(stderr, "ERROR: Cannot open %s for write\n", fn);
      return -1;
    }
    return 0;
  }
  void close() {
    if (this->fp) {
      fclose(this->fp);
      this->fp = NULL;
    }
  }
  int write(const char* s) { return fputs(s, this->fp); }
  int writeLine(const char* s) {
    int ret = fputs(s, this->fp);
    fputc('\n', this->fp);
    return (ret + 1);
  }
  int printf(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    int ret = vfprintf(fp, fmt, args);
    va_end(args);
    return ret;
  }

 private:
  FILE* fp;
};  // end TextFileWriter

class GzipFileWriter : public AbstractFileWriter {
 public:
  GzipFileWriter(const std::string& fn, bool append = false) {
    if (this->open(fn.c_str(), append)) {
      fprintf(stderr, "Cannot create gzip file %s\n", fn.c_str());
    }
  }
  virtual ~GzipFileWriter() {
    this->close();
#ifdef IO_DEBUG
    fprintf(stderr, "GzipFileWriter desc()\n");
#endif
  }
  int open(const char* fn, bool append = false) {
    if (append) fprintf(stderr, "Gzip does not support appending.\n");
    this->fp = gzopen(fn, "wb");
    if (!this->fp) {
      fprintf(stderr, "ERROR: Cannot open %s for write\n", fn);
      return -1;
    }
    return 0;
  }
  void close() {
    if (this->fp) {
      gzclose(this->fp);
      this->fp = NULL;
    }
  }
  int write(const char* s) { return gzputs(this->fp, s); }
  int writeLine(const char* s) {
    int ret = gzputs(this->fp, s);
    gzputc(this->fp, '\n');
    return (ret + 1);
  }

 private:
  gzFile fp;
};  // end GzipFileWriter

class Bzip2FileWriter : public AbstractFileWriter {
 public:
  Bzip2FileWriter(const std::string& fn, bool append = false) : bzp(NULL) {
    if (this->open(fn.c_str(), append)) {
      fprintf(stderr, "Cannot create bzip2 file %s\n", fn.c_str());
    }
  }
  virtual ~Bzip2FileWriter() {
    this->close();
#ifdef IO_DEBUG
    fprintf(stderr, "Bzip2FileWriter desc()\n");
#endif
  }
  int open(const char* fn, bool append = false) {
    if (append) fprintf(stderr, "bzip2 does not support appending.\n");
    this->fp = fopen(fn, "wb");
    if (fp == NULL) return -1;

    this->bzp = BZ2_bzWriteOpen(
        &this->bzerror, this->fp, 9, 0,
        30);  // block size is 9, 0 means silent, 30 means working factor
    if (this->bzerror != BZ_OK) {
      BZ2_bzWriteClose(&bzerror, this->bzp, 0, 0, 0);  // 0: abandon, 0: results
                                                       // of # of bytes for
                                                       // input, 0: results of #
                                                       // of bytes outputted.
      fprintf(stderr, "ERROR: Cannot open %s for write\n", fn);
      return -1;
    }
    return 0;
  }
  void close() {
    BZ2_bzWriteClose(&bzerror, this->bzp, 0, 0, 0);
    if (bzerror != BZ_OK) {
    }
    if (this->fp) fclose(this->fp);

    this->bzp = NULL;
    this->fp = NULL;
  }
  int write(const char* s) {
    int ret = strlen(s);
    BZ2_bzWrite(&this->bzerror, this->bzp, (void*)s, ret);
    if (this->bzerror != BZ_OK) {
      this->close();
      return -1;
    }
    return ret;
  }
  int writeLine(const char* s) {
    int ret = strlen(s);
    BZ2_bzWrite(&this->bzerror, this->bzp, (void*)s, ret);
    if (this->bzerror != BZ_OK) {
      this->close();
      return -1;
    }
    char buf[] = "\n";
    BZ2_bzWrite(&this->bzerror, this->bzp, buf, 1);
    if (this->bzerror != BZ_OK) {
      this->close();
      return -1;
    }
    return (ret + 1);
  }

 private:
  FILE* fp;
  BZFILE* bzp;
  int bzerror;
};  // end Bzip2FileWriter

class BGZipFileWriter : public AbstractFileWriter {
 public:
  BGZipFileWriter(const std::string& fn, bool append = false) {
    if (this->open(fn.c_str())) {
      fprintf(stderr, "Cannot create BGzip file %s\n", fn.c_str());
    }
  }
  virtual ~BGZipFileWriter() {
    this->close();
#ifdef IO_DEBUG
    fprintf(stderr, "BGZipFileWriter desc()\n");
#endif
  }
  /**
   * @param append: ignored
   */
  int open(const char* fn, bool append = false);
  void close();
  int write(const char* s);
  int writeLine(const char* s);

 private:
  BGZF* fp;
};  // end BGZipFileWriter

class StdoutWriter : public AbstractFileWriter {
 public:
  int open(const char* fn, bool append = false) { return 0; }
  void close() {}
  int write(const char* s) { return fputs(s, stdout); }
  int writeLine(const char* s) {
    int ret = fputs(s, stdout);
    putchar('\n');
    return (ret + 1);
  }
};

#define DEFAULT_WRITER_BUFFER 4096
class BufferedFileWriter : public AbstractFileWriter {
 public:
  BufferedFileWriter(AbstractFileWriter* f,
                     int bufLen = DEFAULT_WRITER_BUFFER) {
    this->bufLen = DEFAULT_WRITER_BUFFER;
    this->buf = new char[bufLen + 1];  // last char in the buffer is always '\0'
    // that help to use fputs()
    if (!this->buf) {
      fprintf(stderr, "%s:%d Cannot create BufferedFileWriter\n", __FILE__,
              __LINE__);
      exit(1);
    }
    this->buf[bufLen] = '\0';
    this->bufPtr = 0;
    this->f = f;
  }
  ~BufferedFileWriter() {
    if (this->buf) {
      delete[] this->buf;
      this->buf = NULL;
    }
#ifdef IO_DEBUG
    fprintf(stderr, "BufferedFileWriter desc()\n");
#endif
  }
  int open(const char* fn, bool append = false) {
    return this->f->open(fn, append);
  }
  void close() {
    this->flush();
    // this->f->close();
  }
  int write(const char* s) {
    int nbyte = 0;
    int i = 0;
    while (s[i] != '\0') {
      this->buf[this->bufPtr++] = s[i++];
      nbyte++;
      if (this->bufPtr == this->bufLen) {
        this->f->write(this->buf);
        this->bufPtr = 0;
      }
    }
    return nbyte;
  }
  int writeLine(const char* s) {
    int ret = this->write(s);
    this->write("\n");
    return (ret + 1);
  }
  int flush() {
    this->buf[this->bufPtr] = '\0';
    this->f->write(this->buf);
    this->bufPtr = 0;
    return 0;
  }

 private:
  char* buf;
  int bufLen;
  int bufPtr;
  AbstractFileWriter* f;
};  // end BufferedFileWriter

// static method
AbstractFileReader* AbstractFileReader::open(const char* fileName) {
  AbstractFileReader* fr = NULL;
  if (!fileName || fileName[0] == '\0') {
    fprintf(stderr, "Empty file name.\n");
    return fr;
  }

#ifdef _USE_KNETFILE
  if (strstr(fileName, "ftp://") == fileName ||
      strstr(fileName, "http://") == fileName) {
    fr = new KnetFileReader(fileName);
    // fprintf(stderr, "open knetfile %s\n", fileName);
    return fr;
  }
#endif
  // check fileName suffix
  size_t l = strlen(fileName);
  if (l > 3 && !strcmp(fileName + l - 3, ".gz")) {
    fr = new GzipFileReader(fileName);
    return fr;
  } else if (l > 4 && !strcmp(fileName + l - 4, ".bz2")) {
    fr = new Bzip2FileReader(fileName);
    return fr;
  }

  switch (AbstractFileReader::checkFileType(fileName)) {
    case PLAIN:
      fr = new PlainFileReader(fileName);
      break;
    case GZIP:
      fr = new GzipFileReader(fileName);
      break;
    case BZIP2:
      fr = new Bzip2FileReader(fileName);
      break;
    default:
      fprintf(stderr, "Cannot detect file type (does it exist?!)\n");
      break;
  }
  return fr;
}
// static method
void AbstractFileReader::close(AbstractFileReader** f) {
  assert(f && *f);
  (*f)->close();
  delete (*f);
  *f = NULL;
};

// check header for known file type
FileType AbstractFileReader::checkFileType(const char* fileName) {
  // treat stdin as plain text file
  if (strncmp(fileName, "-", 1) == 0) {
    return PLAIN;
  }
  // magic numbers
  const int gz_magic[2] = {0x1f, 0x8b};  /* gzip magic header */
  const int bzip2_magic[2] = {'B', 'Z'}; /* bzip2 magic header */
  // read file header
  FILE* fp = fopen(fileName, "rb");
  if (!fp) return UNKNOWN;
  unsigned char header[2] = {0, 0};
  int n = fread(header, sizeof(char), 2, fp);
  fclose(fp);
  // check file types
  if (n >= 2 && header[0] == gz_magic[0] && header[1] == gz_magic[1]) {
    return GZIP;
  }
  if (n >= 2 && header[0] == bzip2_magic[0] && header[1] == bzip2_magic[1]) {
    return BZIP2;
  }
  return PLAIN;
  /* // check the characters fall into visible ASCII range */
  /* if ( header[0] >= 0x20 /\* space *\/ && */
  /*      header[0] <  0x7f /\* DEL *\/   && */
  /*      header[1] >= 0x20 /\* space *\/ && */
  /*      header[1] <  0x7f /\* DEL *\/) { */
  /*     return PLAIN; */
  /* }  */
  /* return UNKNOWN; */
};

/**
 * @return number of empty elements filtered out
 */
int removeEmptyField(std::vector<std::string>* fields) {
  int s = fields->size();
  fields->erase(std::remove(fields->begin(), fields->end(), ""), fields->end());
  s -= fields->size();
  return s;
};

AbstractFileWriter::~AbstractFileWriter() {
#ifdef IO_DEBUG
  fprintf(stderr, "AbstractFileWriter desc()\n");
#endif
};

int BGZipFileWriter::open(const char* fn, bool append) {
  if (append) fprintf(stderr, "Gzip does not support appending.\n");
  this->fp = bgzf_open(fn, "w");
  if (!this->fp) {
    fprintf(stderr, "ERROR: Cannot open %s for write\n", fn);
    return -1;
  }
  return 0;
}
void BGZipFileWriter::close() {
  if (this->fp) {
    bgzf_close(this->fp);
    this->fp = NULL;
  }
};
int BGZipFileWriter::write(const char* s) {
  return bgzf_write(this->fp, s, strlen(s));
};
int BGZipFileWriter::writeLine(const char* s) {
  int ret = bgzf_write(this->fp, s, strlen(s));
  ret += bgzf_write(this->fp, "\n", 1);
  return (ret);
};

bool fileExists(std::string fn) {
  FILE* fp = fopen(fn.c_str(), "r");
  if (fp != NULL) {
    fclose(fp);
    return true;
  }

  return false;
}

FileWriter::FileWriter(const std::string& fileName, bool append) {
  if (fileName == "stdout") {
    this->fp = new StdoutWriter();
    this->fpRaw = NULL;
    this->createBuffer();
    return;
  }

  if (this->checkSuffix(fileName.c_str(), ".gz")) {
    this->fpRaw = new GzipFileWriter(fileName.c_str(), append);
  } else if (this->checkSuffix(fileName.c_str(), ".bz2")) {
    this->fpRaw = new Bzip2FileWriter(fileName.c_str(), append);
  } else {
    this->fpRaw = new TextFileWriter(fileName.c_str(), append);
  }
  this->fp = new BufferedFileWriter(this->fpRaw);
  if (!this->fpRaw || !this->fp) {
    fprintf(stderr, "%s:%d Cannot create file\n", __FILE__, __LINE__);
    exit(1);
  }

  this->createBuffer();
}

FileWriter::FileWriter(const std::string& fileName, FileType t) {
  if (fileName == "stdout") {
    this->fp = new StdoutWriter();
    this->fpRaw = NULL;
    this->createBuffer();
    return;
  }

  bool append = false;
  if (PLAIN == t) {
    this->fpRaw = new TextFileWriter(fileName, append);
  } else if (GZIP == t) {
    this->fpRaw = new GzipFileWriter(fileName, append);
  } else if (BZIP2 == t) {
    this->fpRaw = new Bzip2FileWriter(fileName, append);
  } else if (BGZIP == t) {
    this->fpRaw = new BGZipFileWriter(fileName, append);
  } else {
    fprintf(stderr, "Unrecognized file type, use plain text format instead!\n");
    this->fpRaw = new TextFileWriter(fileName, append);
  }

  this->fp = new BufferedFileWriter(this->fpRaw);
  if (!this->fpRaw || !this->fp) {
    fprintf(stderr, "%s:%d Cannot create file\n", __FILE__, __LINE__);
    exit(1);
  }

  this->createBuffer();
}
