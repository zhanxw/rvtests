#include "VCFIndividual.h"

#include "base/IO.h"

char VCFIndividual::defaultValue[2] = {'.', '\0'};
VCFValue VCFIndividual::defaultVCFValue(&(defaultValue[0]), 0, 0);

void VCFIndividual::output(FileWriter* fp) const {
  for (unsigned int i = 0; i < fdLen; ++i) {
    if (i) fp->write(':');
    this->fd[i].output(fp);
  }
}

void VCFIndividual::toStr(std::string* fp) const {
  std::string s;
  for (unsigned int i = 0; i < fdLen; ++i) {
    if (i) fp->push_back(':');
    this->fd[i].toStr(&s);
    fp->append(s);
  }
}
