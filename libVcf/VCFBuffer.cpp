#include "VCFBuffer.h"

#include "base/IO.h"

void VCFBuffer::output(FileWriter* fp, char c) const{
  for (size_t i = 0; i != len; ++i) {
    fp->write(buf[i]);
  }
  fp->write(c);
}
