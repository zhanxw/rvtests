#ifndef _MODELPARSER_H_
#define _MODELPARSER_H_

#include "base/Logger.h"

extern Logger* logger;

/**
 * Parse "mb" to {"mb"}
 * Parse "mb(nperm=10000, alpha=0.1)" then later use assign("nperm", &p) to set value for p 
 */
class ModelParser{
 public:
  /**
   * @return 0: if parse succeed.
   * NOTE: will convert all letters to lower case
   */
  int parse(const char* s){
    std::string arg = s;
    tolower(&arg);

    size_t l = arg.find(ModelParser::LEFT_DELIM);
    if ( l == std::string::npos){
      this->name = arg;
      return 0;
    }
    this->name = arg.substr(0, l);
    if (arg[arg.size() - 1] != ModelParser::RIGHT_DELIM){
      logger->error("Please use this format: model(model_param1=v1)");
      return -1;
    }
    std::vector<std::string> params;
    std::string allParam = arg.substr(l + 1, arg.size() - 1 - 1 -l);
    int ret = stringTokenize(allParam, ',', &params);
    UNUSED(ret);
    for (size_t i = 0; i < params.size(); ++i) {
      l = params[i].find('=');
      if (l == std::string::npos) {
        this->param[params[i]] = "";
      } else {
        std::string key = params[i].substr(0, l);
        std::string value = params[i].substr(l + 1, params[i].size() - l - 1);
        this->param[key] = value;
      };
    };
    logger->info("load %zu parameters.", this->size());
    return 0;
  }
  int parse(std::string& s){
    return this->parse(s.c_str());
  }
  const std::string& getName() const {
    return this->name;
  };
  /**
   * assign @param tag to @param value
   */
  ModelParser& assign(const char* tag, bool* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = true;
    } else {
      (*value) = false;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, double* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atof(this->param[tag]);
    } else {
      logger->error("Cannot find parameter [ %s ]", tag);
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, int* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atoi(this->param[tag]);
    } else {
      logger->error("Cannot find parameter [ %s ]", tag);
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, std::string* value){
    if (this->param.find(tag) != this->param.end()){
      (*value) = (this->param[tag]);
    } else {
      logger->error("Cannot find parameter [ %s ]", tag);
    }
    return (*this);
  };
  /**
   * if @param tag is provided, assign @param tag to @param value
   * otherwise, set @param value to @param def
   */
  ModelParser& assign(const char* tag, bool* value, const bool def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = true;
    } else {
      (*value) = def;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, double* value, const double def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atof(this->param[tag]);
    } else {
      (*value) = def;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, int* value, const int def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = atoi(this->param[tag]);
    } else {
      (*value) = def;
    }
    return (*this);
  };
  ModelParser& assign(const char* tag, std::string* value, const std::string& def){
    if (this->param.find(tag) != this->param.end()){
      (*value) = (this->param[tag]);
    } else {
      (*value) = def;
    }
    return (*this);
  };
  const size_t size() const {
    return this->param.size();
  };
  static const char LEFT_DELIM = '(';
  static const char RIGHT_DELIM = ')';
 private:
  std::string name;
  std::map<std::string, std::string> param;
};



#endif /* _MODELPARSER_H_ */
