#include "VCFValue.h"

#include "base/IO.h"

char VCFValue::defaultValue[1] = {'\0'};

void VCFValue::output(FileWriter* fp) const {
  if (!line) return;
  for (int i = beg; i < end; ++i) {
    fp->write(line[i]);
  }
};

void VCFValue::output(FileWriter* fp, char c) const {
  for (int i = beg; i <= end; ++i) {
    fp->write(line[i]);
  }
  fp->write(c);
}
