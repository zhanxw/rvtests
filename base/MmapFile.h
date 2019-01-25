#ifndef _MMAPFILE_H_
#define _MMAPFILE_H_

#ifdef _WIN32
#include <io.h>
#include <stdint.h>
#include <stdio.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#endif

extern size_t getFileSize(const char* fileName);

class MmapFile {
 public:
  MmapFile();
  MmapFile(const char* fileName);
  virtual ~MmapFile();

  int open(const char* fileName);
  int close();
  void* data;  // mmap() data goes here
  size_t getFileSize() { return this->fileSize; };

 private:
  size_t fileSize;
#ifdef _WIN32
  HANDLE handle;
#endif
};

#endif /* _MMAPFILE_H_ */
