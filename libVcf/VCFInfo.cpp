/#include "VCFInfo.h"

const VCFValue VCFInfo::defaultValue;

const VCFValue& VCFInfo::getTag(const char* tag, bool* exists){
  if (!tag || tag[0] == '\0') {
    *exists = false;
    return VCFInfo::defaultValue;
  }

  if (!this->hasParsed)
    this->parseActual();

  std::string s = tag;
  if (!this->data.find(s) ) {
    *exists = false;
    return VCFInfo::defaultValue;
  } else {
    *exists = true;
    return this->data[s];
  };
};
