#include "MmapFile.h"

#include <fcntl.h>  // O_RDONLY
#define __STDC_LIMIT_MACROS
#define __STDC_CONSTANT_MACROS
#include <stdint.h>  // SIZE_MAX

#ifdef _WIN32
#include <io.h>
#include <stdio.h>
#include <sys/stat.h>
#include <windows.h>
#else
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#endif

MmapFile::MmapFile() : data(0) {}

MmapFile::MmapFile(const char* fileName) : data(0) { open(fileName); }

MmapFile::~MmapFile() { close(); }

size_t getFileSize(const char* fileName) {
  struct stat fileStat;
  int fd = open(fileName, O_RDONLY);
  if (fd == -1) {
    fprintf(stderr, "Cannot open file");
    // exit(1);
    return SIZE_MAX;
  }

  if (fstat(fd, &fileStat)) {
    fprintf(stderr, "Cannot fstat() file");
    // exit(1);
    return SIZE_MAX;
  }
  close(fd);
  return fileStat.st_size;
}

int MmapFile::open(const char* fileName) {
  int filedes = ::open(fileName, O_RDONLY);
  if (filedes < 0) {
    fprintf(stderr, "Cannot open file");
    // exit(1);
    return -1;
  }
  this->fileSize = ::getFileSize(fileName);
  if (data) {
    this->close();
  }
#ifndef _WIN32
  this->data = mmap(0, this->fileSize, PROT_READ, MAP_SHARED, filedes, 0);
  if (this->data == MAP_FAILED) {
    fprintf(stderr, "mmap() failed!");
    // exit(1);
    return -1;
  }
#else
  // see:
  // http://courses.washington.edu/hypertxt/cwb/cl/windows-mmap.c
  // https://docs.microsoft.com/en-us/windows/desktop/api/winbase/nf-winbase-createfilemappinga
  this->handle = CreateFileMapping(
      (HANDLE)_get_osfhandle(filedes),
      NULL,  //  the file mapping object gets a default security descriptor
      PAGE_READONLY,  // can also be PAGE_WRITECOPY,
      0,  // dwMaximumSizeHigh are 0 (zero), the maximum size of the file
          // mapping object is equal to the current size of the file that hFile
          // identifies
      0,  // dwMaximumSizeLow (see above)
      NULL);  // lpName = NULL, the file mapping object is created without a
              // name.
  if (handle == NULL) {
    return -1;
  }
  this->data = MapViewOfFileEx(handle, FILE_MAP_READ,
                               0,  // dwFileOffsetHigh: The high-order DWORD of
                                   // the file offset where the view is to begin
                               0,  // dwFileOffsetLow: see above
                               this->fileSize,  // dwNumberOfBytesToMap
                               NULL);  // lpBaseAddress: A pointer to the memory
                                       // address in the calling process address
                                       // space where mapping begins.

#endif
  return 0;
}

int MmapFile::close() {
#ifndef _WIN32
  // @return -1: means error
  if (munmap(this->data, this->fileSize)) {
    return -1;
  }
#else
  if (!UnmapViewOfFile(this->data)) {
    return -1;
  }
  if (!CloseHandle(this->handle)) {
    fprintf(stderr, "unable to close file mapping handle\n");
    return -1;
  }

#endif
  this->data = 0;
  return 0;
}
