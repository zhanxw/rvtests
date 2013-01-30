#ifndef _RESULT_H_
#define _RESULT_H_

#include "base/TypeConversion.h"

/**
 * Store key-value pair, 
 */
class Result{
public:
  Result () : defaultValue("NA") {};
  void addHeader(const char* key) {
    data[key] = defaultValue;
  }
  void addHeader(const std::string& key) {
    data[key] = defaultValue;
  }

  bool existHeader(const char* key) {
    if (!data.find(key)) {
      fprintf(stderr, "Cannot find [ %s ] in result header...\n", key);
      return false;
    }
    return true;
  }
  bool existHeader(const std::string& key) {
    return existHeader(key.c_str());
  }
  void updateValue(const char* key, const char* val) {
    if (!existHeader(key)) {
      return;
    }
    data[key] = val;
  }
  void updateValue(const char* key, const std::string& val) {
    if (!existHeader(key)) {
      return;
    }
    data[key] = val;    
  }
  void updateValue(const std::string& key, const std::string& val) {
    if (!existHeader(key)) {
      return;
    }
    data[key] = val;
  }
  void updateValue(const char* key, const int val) {
    if (!existHeader(key)) {
      return;
    }
    data[key] = toString(val);
  }
  void updateValue(const char* key, const double val) {
    if (!existHeader(key)) {
      return;
    }
    data[key] = floatToString(val);
  }
  
  void clearValue() {
    int n = data.size();
    for (int i = 0; i < n ; ++i ) {
      this->data.valueAt(i) = defaultValue;
    }
  }
  /**
   * Write the keys separated by '\t'
   */
  void writeHeader(FILE* fp) const {
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      if (i){
        fputc('\t', fp);
      }
      fputs(data.keyAt(i).c_str(), fp);
    }
  }
  void writeHeaderTab(FILE* fp) const {
    writeHeader(fp);
    fputc('\t', fp);    
  }

  void writeHeaderLine(FILE* fp) const {
    writeHeader(fp);
    fputc('\n', fp);    
  }
  
  /**
   * Write the values separated by '\t'
   */
  void writeValue(FILE* fp) const{
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      if (i){
        fputc('\t', fp);
      }
      fputs(data.valueAt(i).c_str(), fp);
    }
  }
  void writeValueTab(FILE* fp) const{
    writeValue(fp);
    fputc('\t', fp);    
  }
  void writeValueLine(FILE* fp) const{
    writeValue(fp);
    fputc('\n', fp);    
  }
  /**
   * Use '\t' to join headers
   */
  std::string joinHeader() const {
    std::string s;
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      s += data.keyAt(i);
      if (i){
        s += '\t';
      }
    }
    return s;
  }

  std::string joinValue() const{
    std::string s;
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      s += data.valueAt(i);
      if (i){
        s += '\t';
      }
    }
    return s;
  }
  
  const std::string& operator[] (const std::string& key) const{
    if (data.find(key)) {
      return data[key];
    }
    return defaultValue;
  }
  const std::string& operator[] (const char* key) const{
    if (data.find(key)) {
      return data[key];
    }
    return defaultValue;
  }
private:
  OrderedMap<std::string, std::string> data;
  std::string defaultValue;
};

#endif /* _RESULT_H_ */
