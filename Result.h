#ifndef _RESULT_H_
#define _RESULT_H_

#include "base/OrderedMap.h"
#include "base/TypeConversion.h"
#include "base/IO.h"

/**
 * Store key-value pair for minimal typing
 * Internally, all keys and values are strings
 * Memory layout:
 *  key1 -> val, val, val....
 *  key2 -> val, val, val....
 *           ^    ^    ^
 *           |    |    |
 *           L1   L2   L3
 *           We assume updating data are performed layer by layer
 *           If cross layout updating value happened, we will generate an error.
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

  //////////////////////////////////////////////////
  // Use FILE* to output
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
//////////////////////////////////////////////////
  // Use FileWriter* to output
  /**
   * Write the keys separated by '\t'
   */
  void writeHeader(FileWriter* fp) const {
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      if (i){
        fp->write('\t');
      }
      fp->write(data.keyAt(i).c_str());
    }
  }
  void writeHeaderTab(FileWriter* fp) const {
    writeHeader(fp);
    fp->write('\t');
  }

  void writeHeaderLine(FileWriter* fp) const {
    writeHeader(fp);
    fp->write('\n');    
  }
  
  /**
   * Write the values separated by '\t'
   */
  void writeValue(FileWriter* fp) const{
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      if (i){
        fp->write('\t');
      }
      fp->write(data.valueAt(i).c_str());
    }
  }
  void writeValueTab(FileWriter* fp) const{
    writeValue(fp);
    fp->write('\t');
  }
  void writeValueLine(FileWriter* fp) const{
    writeValue(fp);
    fp->write('\n');
  }
  
  /**
   * Use '\t' to join headers
   */
  std::string joinHeader() const {
    std::string s;
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      if (i){
        s += '\t';
      }
      s += data.keyAt(i);
    }
    return s;
  }

  std::string joinValue(const char c = '\t') const{
    std::string s;
    int n = data.size();
    for (int i = 0; i < n; ++i) {
      if (i){
        s += c;
      }
      s += data.valueAt(i);
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
