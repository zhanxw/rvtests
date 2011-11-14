#ifndef _IO_H_
#define _IO_H_

#include <stdio.h>  //fopen
#include <stdlib.h> //malloc
#include <string.h> //strchr
#include <assert.h> //assert
#include <stdarg.h> // va_list

#include <string>
#include <vector>

// #define IO_DEBUG

/**
 * Sample usage:
 * std::string line;
 * FileReader* f = FileReader::open("abc");
 * while (f->readLine(line) > 0) {
 *   ...
 * }
 * FileReader::close(&f);
 */
class FileReader{
  public:
    typedef enum FileType {
        PLAIN = 0,
        GZIP = 1,
        BZIP2 = 2,
        UNKNOWN = 99
    } FileType;
    virtual ~FileReader() {}; // make it virtual so subclass types can close file handle
    static FileReader* open(const char* fileName);
    static void close(FileReader** f);
    // virtual functions
    // each specific file type will need to implement the following function
    //virtual int readLine(std::string* line) = 0;
    //virtual int readLineBySep(std::vector<std::string>* fields, const char* sep) = 0;
    virtual int getc() = 0;
    virtual bool isEof() = 0;
    virtual void close() = 0;
    virtual int read(void* buf, unsigned int len) = 0;
    // common utility function
    static FileType checkFileType(const char* fileName);
  protected:
    FileReader() {}; // forbid explicit create FileReader class.
};

class PlainFileReader: public FileReader{
  public:
  PlainFileReader(const char* fileName):
    fp(NULL) {
        this->open(fileName);
#ifdef IO_DEBUG
        fprintf(stderr, "PlainFileReader() open %s\n", fileName);
#endif
    };
    virtual ~PlainFileReader() {
#ifdef IO_DEBUG
        fprintf(stderr, "~PlainFileReader() close\n");
#endif
        this->close();
    };
    
    // get a char, if EOF, return EOF
    int getc(){
        return ::getc(this->fp);
    }
    // check eof 
    bool isEof() {
        return (feof(this->fp) != 0);
    }
    // open
    FILE* open(const char* fileName) {
        this->fp = fopen(fileName, "r");
        if (!this->fp) {
            fprintf(stderr, "ERROR: Cannot open %s\n", fileName);
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
    int read(void* buf, unsigned int len) {
        return ::fread(buf, sizeof(char), len, this->fp);
    };
  private:
    FILE* fp;
};

//////////////////////////////////////////////////////////////////////
// Gzip reading class
#include <zlib.h>
class GzipFileReader: public FileReader{
  public:
  GzipFileReader(const char* fileName):
    fp(NULL) {
        this->open(fileName);
#ifdef IO_DEBUG
        fprintf(stderr, "GzipFileReader() open %s\n", fileName);
#endif
    };
    virtual ~GzipFileReader() {
#ifdef IO_DEBUG
        fprintf(stderr, "~PlainFileReader() close\n");
#endif
        this->close();
    };
    
    // get a char, if EOF, return EOF
    int getc(){
        return gzgetc(this->fp);
    }
    // check eof 
    bool isEof() {
        return (gzeof(this->fp) != 0);
    }
    // open
    gzFile open(const char* fileName) {
        this->fp = gzopen(fileName, "r");
        if (!this->fp) {
            fprintf(stderr, "ERROR: Cannot open %s\n", fileName);
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
    int read(void* buf, unsigned int len) {
        return gzread(this->fp, buf, len);
    };

  private:
    gzFile fp;
};
//////////////////////////////////////////////////////////////////////
// Bzip2 reading class
#include <bzlib.h>
class Bzip2FileReader: public FileReader{
  public:
  Bzip2FileReader(const char* fileName):
    fp(NULL) {
        this->open(fileName);
#ifdef IO_DEBUG
        fprintf(stderr, "Bzip2FileReader() open %s\n", fileName);
#endif
    };
    virtual ~Bzip2FileReader() {
#ifdef IO_DEBUG
        fprintf(stderr, "~Bzip2FileReader() close\n");
#endif
        if (this->fp) {
            BZ2_bzclose(fp);
        }
    };
    
    // get a char, if EOF, return EOF
    int getc(){
        char c;
        this->bzerror = BZ_OK;
        int nBuf = BZ2_bzRead(&this->bzerror, this->bzp, &c, sizeof(char));
        if (this->bzerror == BZ_OK) {
            return c;
        } else {
            return EOF;
        }
    }
    // check eof 
    bool isEof() {
        return (this->bzerror == BZ_STREAM_END);
    }
    // open
    BZFILE* open(const char* fileName) {
        this->fp = fopen(fileName, "r");
        if (!this->fp) {
            fprintf(stderr, "ERROR: Cannot open %s\n", fileName);
            return NULL;
        }
        this->bzp = BZ2_bzReadOpen(&this->bzerror, this->fp, 0, 0, NULL, 0);

        if (this->bzerror != BZ_OK) {
            BZ2_bzReadClose ( &bzerror, this->bzp );
            fprintf(stderr, "ERROR: Cannot open %s\n", fileName);
            return NULL;
        }
    }
    // close 
    void close() {
        if (this->bzerror != BZ_STREAM_END) {
            BZ2_bzReadClose(&this->bzerror, fp);
            /* omit code to handle error */
        } else {
            BZ2_bzReadClose(&this->bzerror, fp);
        }
        this->fp = NULL;
        this->bzp = NULL;
        this->bzerror = 0;
    };
    int read(void* buf, unsigned int len) {
        return BZ2_bzRead ( &this->bzerror, this->bzp, buf, len);
    };
  private:
    FILE* fp;
    BZFILE* bzp;
    int bzerror;
};
//////////////////////////////////////////////////////////////////////

//static method
FileReader* FileReader::open(const char* fileName){
    FileReader * fr = NULL;
    switch(FileReader::checkFileType(fileName)) {
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
        fprintf(stderr, "Cannot detect file type (even not text file?!\n");
        break;
    }
    return fr;
}
// static method
void FileReader::close(FileReader** f) {
    assert(f && *f);
    (*f)->close();
    delete (*f);
    *f = NULL;
};

// check header for known file type
FileReader::FileType FileReader::checkFileType(const char* fileName){
    // treat stdin as plain text file
    if (strncmp(fileName, "-", 1) == 0) {
        return PLAIN;
    }
    // magic numbers
    const int gz_magic[2] = {0x1f, 0x8b}; /* gzip magic header */
    const int bzip2_magic[2] = {'B', 'Z'}; /* bzip2 magic header */
    // read file header    
    FILE* fp = fopen(fileName, "rb");
    if (!fp) return UNKNOWN;
    unsigned char header[2]={0,0};
    int n = fread(header, sizeof(char), 2, fp);
    fclose(fp);
    // check file types
    if ( n >= 2 && header[0] == gz_magic[0] && header[1] == gz_magic[1]) {
        return GZIP;
    }
    if ( n >= 2 && header[0] == bzip2_magic[0] && header[1] == bzip2_magic[1]) {
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
class BufferedReader: public FileReader{
  public:
    BufferedReader(const char* fileName, int bufferCapacity) {
#ifdef IO_DEBUG
        fprintf(stderr, "BufferedReader open %s\n", fileName);
#endif
        // initialize buf
        if (bufferCapacity <= 0) {
            fprintf(stderr, "Buffer size should be greater than 0, now use default buffer size 8096 instead of %d.\n", bufferCapacity);
            this->bufCap = 8096;
        } else {
            this->bufCap = bufferCapacity;
        }
        this->buf = (char*) malloc( sizeof(char)* this->bufCap);
        if (!this->buf) {
            fprintf(stderr, "Cannot allocate buffer for BufferedReader.\n");
            return;
        }
        this->bufPtr = 0;
        this->bufEnd = 0;
        // initialize fp
        this->fp = FileReader::open(fileName);
        if (!this->fp) {
            fprintf(stderr, "Canont open file %s\n", fileName);
            this->fp = NULL;
        }
    }
    virtual ~BufferedReader(){
        this->close();
    }
    int getc() {
        if (this->bufPtr == this->bufEnd) { //buffer all used, need to refresh
            this->bufEnd = this->fp->read(this->buf, this->bufCap);
            this->bufPtr = 0;
        } 

        if (this->bufPtr < this->bufEnd)
            return (this->buf [ this->bufPtr++] );
        else 
            return EOF;
    }
    bool isEof() {
        if (this->fp->isEof() && this->bufPtr == this->bufEnd){
            return true;
        }
        return false;
    }
    void close() {
#ifdef IO_DEBUG
        fprintf(stderr, "BufferedReader close\n");
#endif
        // delete fp
        if (this->fp){
            FileReader::close(&fp);
        }
        // delete buf
        if (this->buf) {
            free(this->buf);
            this->buf = 0;
            this->bufCap = 0;
            this->bufPtr = 0;
            this->bufEnd = 0;
        }
    }
    int read(void* buf, unsigned int len) {
        // use current buffer to fill in buf
        unsigned int idx = 0;
        while (this->bufPtr < this->bufEnd && len > 0) {
            ((char*)buf)[idx++] = this->buf[this->bufPtr++];
            len --;
        }
        if (len == 0) {
            return idx;
        }
        // fill rest of buf
        unsigned int nRead = this->fp->read(((char*)buf)+idx, len);
        idx += nRead;
        // refill buffer
        this->bufEnd = this->fp->read(this->buf, this->bufCap);
        this->bufPtr = 0;
        return idx;
    }
  private:
    unsigned int bufCap; // capacity of the buffer
    unsigned int bufEnd; // bufPtr should not read beyond bufEnd(incluive)
    unsigned int bufPtr; // from which buffer begins to read
    char* buf;
    FileReader* fp;
};

/** Example code:
// LineReader lr(fn);
// while(lr.readLine(&line)>0){
//         fprintf(stdout, "%s\n", line.c_str());
// }
*/
class LineReader{
  public:
    LineReader(const char* fileName){
#ifdef IO_DEBUG
        fprintf(stderr, "LineReader open %s\n", fileName);
#endif
        this->fp = new BufferedReader(fileName, 1024);
        if (!this->fp) {
            fprintf(stderr, "Canont open file %s\n", fileName);
            this->fp = NULL;
        }
    }
    LineReader(FileReader* fp) {
        this->fp = fp;    
    }
    virtual ~LineReader(){
        if (this->fp){
            fp->close();
            delete this->fp;
            this->fp = NULL;
        }
#ifdef IO_DEBUG
        fprintf(stderr, "LineReader close\n");
#endif
    }
    // return number of characters read.
    // when reading an empty line, will return 1, as we read '\n', however, line will be empty
    // when reading the end, we will return 0
    int readLine(std::string* line) {
        assert(this->fp && line);
        if (this->fp->isEof()) return 0;
        line->clear();
        char c;
        unsigned nRead = 0;
        while (true) {
            c = this->fp->getc();
            if (c == EOF) {
                return nRead;
            } else if (c == '\r') {
                // skip this
                continue;
            } else if (c == '\n') {
                ++nRead;
                return nRead;
            } else { // normal characters
                ++nRead;
                line->push_back(c);
            }
        }   
        assert(false); // should not reach here
        return 0;
    };
    // return number of fields read.
    // when reading an empty line, will return 1, meaning 1 field are read, although its content is empty
    // consecutive separators, e.g. \t\t, will yield empty field
    // when reading to the EOF, will 0. 
    int readLineBySep(std::vector<std::string>* fields, const char* sep) {
        assert(this->fp && fields && sep);
        if (this->fp->isEof()) return 0;
        fields->clear();
        char c;
        std::string s;
        while (true) {
            c = this->fp->getc();
            if (c == EOF) {
                fields->push_back(s);
                return fields->size();
            } else if (c == '\r') {
                // skip this
                continue;
            } else if (c == '\n') {
                fields->push_back(s);
                return fields->size();
            } else if (strchr(sep, c) != NULL) { // separator
                fields->push_back(s);
                s.clear();
            } else { // normal characters
                s.push_back(c);
            }
        }   
        assert(false); // should not reach here
        return 0;
    };
  private:
    FileReader* fp;
};

// FileWriter related classes
class AbstractFileWriter{
  public:
    /// when open is successful, return 0; else: return non-zero
    virtual int open(const char* fn, bool append = false) = 0;
    virtual void close() = 0;
    virtual int write(const char* s) = 0;
    virtual int writeLine(const char* s) = 0;
    virtual ~AbstractFileWriter() = 0;
};
AbstractFileWriter::~AbstractFileWriter() {
#ifdef IO_DEBUG
    fprintf(stderr, "AbstractFileWriter desc()\n");
#endif
};

class TextFileWriter:public AbstractFileWriter{
  public:    
    TextFileWriter(const char* fn, bool append = false){
        if (this->open(fn, append)){
            fprintf(stderr, "Cannot create text file %s\n", fn);
        }
    }
    virtual ~TextFileWriter(){
#ifdef IO_DEBUG
        fprintf(stderr, "TextFileWriter desc()\n");
#endif
        this->close();
    }
    int open(const char* fn, bool append = false){
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
    void close(){
        if (this->fp) {
            fclose(this->fp);
            this->fp = NULL;
        }
    };
    int write(const char* s) {
        return fputs(s, this->fp);
    };
    int writeLine(const char* s) {
        int ret = fputs(s, this->fp);
        fputc('\n', this->fp);
        return (ret + 1);
    };
    int printf(const char *fmt, ...){
        va_list args;
        va_start(args,fmt);
        int ret = vfprintf(fp, fmt, args);
        va_end(args);
        return ret;
    };
  private:
    FILE* fp;
}; // end TextFileWriter

class GzipFileWriter:public AbstractFileWriter{
  public:    
    GzipFileWriter(const char* fn, bool append = false){
        if (this->open(fn, append)){
            fprintf(stderr, "Cannot create text file %s\n", fn);
        }
    }
    virtual ~GzipFileWriter(){
        this->close();
#ifdef IO_DEBUG
        fprintf(stderr, "GzipFileWriter desc()\n");
#endif
    };
    int open(const char* fn, bool append = false){
        if (append) 
            fprintf(stderr, "Gzip does not support appending.\n");
        this->fp = gzopen(fn, "wb");
        if (!this->fp) {
            fprintf(stderr, "ERROR: Cannot open %s for write\n", fn);
            return -1;
        }
        return 0;
    }
    void close(){
        if (this->fp) gzclose(fp);
    };
    int write(const char* s) {
        return gzputs(this->fp, s);
    };
    int writeLine(const char* s) {
        int ret = gzputs(this->fp, s);
        gzputc(this->fp, '\n');
        return (ret + 1);
    };
  private:
    gzFile fp;
}; // end GzipFileWriter

#define DEFAULT_WRITER_BUFFER 4096
class BufferedFileWriter: public AbstractFileWriter{
  public:
    BufferedFileWriter(AbstractFileWriter* f, unsigned int bufLen = DEFAULT_WRITER_BUFFER){
        this->bufLen = DEFAULT_WRITER_BUFFER;
        this->buf = new char[bufLen + 1]; // last char in the buffer is always '\0'
                                          // that help to use fputs()
        if (!this->buf) {
            fprintf(stderr, "Cannot create BufferedFileWriter\n");
            abort();
        }
        this->buf[bufLen] = '\0';
        this->bufPtr = 0;

        if (!this->buf) {
            fprintf(stderr, "Buffer allocation failed!\n");
        }
        this->f = f;
    }
    ~BufferedFileWriter(){
        if (this->buf) {
            delete [] this->buf;
            this->buf = NULL;
        }
#ifdef IO_DEBUG
        fprintf(stderr, "BufferedFileWriter desc()\n");
#endif
    };
    int open(const char* fn, bool append = false){
        return this->f->open(fn, append);
    };
    void close() {
        this->flush();
        this->f->close();
    };
    int write(const char* s){
        int nbyte = 0;
        int i = 0;
        while (s[i] != '\0'){
            this->buf[this->bufPtr++] = s[i++];
            nbyte ++ ;
            if (this->bufPtr == this->bufLen) {
                this->f->write(this->buf);
                this->bufPtr = 0;
            }
        }
        return nbyte;
    };
    int writeLine(const char* s){
        int ret = this->write(s);
        this->write("\n");
        return (ret + 1) ; 
    };
    int flush() {
        this->buf[this->bufPtr] = '\0';
        this->f->write(this->buf);
        this->bufPtr = 0;
        return 0;
    };
  private:
    char* buf;
    unsigned int bufLen;
    unsigned int bufPtr;
    AbstractFileWriter* f;
}; // end BufferedFileWriter


/**
 * usage:
 * FileWriter* fout = new FileWriter("a.txt", "w");
 * fout->write("abc");
 * fout->writeLn("abc");
 * fout->close();
 * delete fout->write;
 */
class FileWriter{
  public:
    FileWriter(const char* fileName, bool append = false){
        int l = strlen(fileName);
        if (this->checkSuffix(fileName, ".gz")) {
            this->fpRaw = new GzipFileWriter(fileName, append);
        } else {
            this->fpRaw = new TextFileWriter(fileName, append);
        }
        this->fp = new BufferedFileWriter(this->fpRaw);
        if (!this->fpRaw || !this->fp){
            fprintf(stderr, "Cannot create file\n");
            abort();
        }

        // create buffer for formatted string
        this->bufLen = 1024;
        this->buf = new char[this->bufLen];
        if (!this->buf){
            fprintf(stderr, "Cannot allocate printf buffer for FileWriter.\n");
        };
    }
    ~FileWriter(){
        if (this->fp){
            this->fp->close();
            delete this->fp;
            this->fp = NULL;
        }
        if (this->fpRaw){
            delete this->fpRaw;
            this->fpRaw = NULL;
        }
        if (this->buf){
            delete [] this->buf;
            this->buf = NULL;
        }
#ifdef IO_DEBUG
        fprintf(stderr, "FileWriter desc()\n");
#endif
    };
    int write(const char* s){
        return this->fp->write(s);
    };
    int writeLine(const char* s){
        int ret = this->fp->write(s);
        this->fp->write("\n");
        return (ret + 1);
    };
    // if @param fileName ends with @param suffix, then return true;
    bool checkSuffix(const char* fileName, const char* suffix){
        int lf = strlen(fileName);
        int ls = strlen(suffix);
        if (lf < ls) return false;
        for (int i = lf - ls, j = 0; j < ls;){
            if (fileName[i++] != suffix[j++]) return false;
        }
        return true;
    };

    int printf(const char *fmt, ...){
        // we'll put the formatted string  to internal buffer
        va_list args;
        va_start(args,fmt);
        int ret;
        while (( ret = vsnprintf(this->buf, this->bufLen, fmt, args)) < 0){
            this->doubleBuffer();
        };
        va_end(args);
        return this->write(this->buf);
    };

  private:
    void doubleBuffer(){
        delete [] this->buf;
        this->buf = NULL;

        this->bufLen *= 2;
        this->buf = new char[bufLen];
        if (!this->buf) {
            fprintf(stderr, "Cannot reallocate printf buffer for FileWriter.\n");
            abort();
        }
    };
    AbstractFileWriter* fp;
    AbstractFileWriter* fpRaw;
    char* buf;
    unsigned int bufLen;
};

#endif /* _IO_H_ */

