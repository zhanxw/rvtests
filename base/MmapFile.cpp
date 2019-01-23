#include "MmapFile.h"

size_t getFileSize(const char* fileName) {
  struct stat fileStat;
  int fd = open(fileName, O_RDONLY);
  if (fd == -1) {
    perror("Cannot open file");
    exit(1);
  }

  if (fstat(fd, &fileStat)) {
    perror("Cannot fstat() file");
    exit(1);
  }
  close(fd);
  return fileStat.st_size;
};

int MmapFile::open(const char* fileName) {
  this->filedes = ::open(fileName, O_RDONLY);
  if (filedes < 0) {
    perror("Cannot open file");
    exit(1);
  }
  this->fileSize = ::getFileSize(fileName);
  if (data) {
    this->close();
  }
  this->data = mmap(0, this->fileSize, PROT_READ, MAP_SHARED, this->filedes, 0);
  if (this->data == MAP_FAILED) {
    perror("mmap() failed!");
    exit(1);
  }
  return 0;
};
void MmapFile::close() {
  munmap(this->data, this->fileSize);
  this->data = 0;
};
