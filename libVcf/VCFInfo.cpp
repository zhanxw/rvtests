#include "VCFInfo.h"

const VCFValue VCFInfo::defaultValue;

const VCFValue& VCFInfo::getTag(const char* tag, bool* isMissing){
  if (!tag || tag[0] == '\0') {
    *isMissing = true;
    return VCFInfo::defaultValue;
  }

  if (!this->hasParsed)
    this->parseActual();

  std::string s = tag;
  if (!this->data.find(s) ) {
    *isMissing = true;
    return VCFInfo::defaultValue;
  } else {
    *isMissing = false;
    return this->data[s];
  };
};
