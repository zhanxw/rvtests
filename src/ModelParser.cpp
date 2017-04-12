#include "ModelParser.h"

#include "base/Logger.h"
#include "base/Utils.h"

#define UNUSED(x) ((void)(x))

extern Logger* logger;

const std::string ModelParser::PARAM_DELIM = ":,";

int ModelParser::parse(const std::string& s) {
  std::string arg = s;
  tolower(&arg);

  size_t l = arg.find(ModelParser::LEFT_DELIM);
  if (l == std::string::npos) {
    this->name = arg;
    return 0;
  }
  this->name = arg.substr(0, l);
  if (arg[arg.size() - 1] != ModelParser::RIGHT_DELIM) {
    logger->error("Please use this format: model(model_param1=v1)");
    return -1;
  }
  std::vector<std::string> params;
  std::string allParam = arg.substr(l + 1, arg.size() - 1 - 1 - l);
  int ret = stringTokenize(allParam, ModelParser::PARAM_DELIM, &params);
  UNUSED(ret);
  for (size_t i = 0; i < params.size(); ++i) {
    l = params[i].find('=');
    if (l == std::string::npos) {
      this->param[params[i]] = "";
    } else {
      std::string key = params[i].substr(0, l);
      std::string value = params[i].substr(l + 1, params[i].size() - l - 1);
      this->param[key] = value;
    }
  }
  logger->info("load %zu parameters.", this->size());
  return 0;
}

const std::string& ModelParser::getName() const { return this->name; }
bool ModelParser::hasTag(const std::string& tag) const {
  return (this->param.find(tolower(tag)) != this->param.end());
}
const char* ModelParser::value(const std::string& tag) const {
  if (hasTag(tag)) {
    return this->param.find(tolower(tag))->second.c_str();
  }
  return NULL;
}
size_t ModelParser::size() const { return this->param.size(); }

/**
 * assign @param tag to @param value
 */
const ModelParser& ModelParser::assign(const std::string& tag,
                                       bool* value) const {
  if (this->hasTag(tag)) {
    (*value) = true;
  } else {
    (*value) = false;
  }
  return (*this);
}
const ModelParser& ModelParser::assign(const std::string& tag,
                                       double* value) const {
  if (this->hasTag(tag)) {
    (*value) = atof(this->value(tag));
  } else {
    logger->error("Cannot find parameter [ %s ]", tag.c_str());
  }
  return (*this);
}
const ModelParser& ModelParser::assign(const std::string& tag,
                                       int* value) const {
  if (this->hasTag(tag)) {
    double d = atof(this->value(tag));
    *value = (int)d;
    if (!isInteger(d)) {
      logger->warn(
          "Convert parameter [ %s ] from specified value [ %s ] to integer [ "
          "%d ]",
          tag.c_str(), this->value(tag), *value);
    }
  } else {
    logger->error("Cannot find parameter [ %s ]", tag.c_str());
  }
  return (*this);
}
const ModelParser& ModelParser::assign(const std::string& tag,
                                       std::string* value) const {
  if (this->hasTag(tag)) {
    (*value) = (this->value(tag));
  } else {
    logger->error("Cannot find parameter [ %s ]", tag.c_str());
  }
  return (*this);
}
/**
 * if @param tag is provided, assign @param tag to @param value
 * otherwise, set @param value to @param def
 */
const ModelParser& ModelParser::assign(const std::string& tag, bool* value,
                                       const bool def) const {
  if (this->hasTag(tag)) {
    (*value) = true;
  } else {
    (*value) = def;
  }
  return (*this);
}
const ModelParser& ModelParser::assign(const std::string& tag, double* value,
                                       const double def) const {
  if (this->hasTag(tag)) {
    (*value) = atof(this->value(tag));
  } else {
    (*value) = def;
  }
  return (*this);
}
const ModelParser& ModelParser::assign(const std::string& tag, int* value,
                                       const int def) const {
  if (this->hasTag(tag)) {
    // convert double then to integer
    // so that scientific notation of float numbers are supported
    double d = atof(this->value(tag));
    *value = (int)d;
    if (!isInteger(d)) {
      logger->warn(
          "Convert parameter [ %s ] from specified value [ %s ] to integer [ "
          "%d ]",
          tag.c_str(), this->value(tag), *value);
    }
  } else {
    (*value) = def;
  }
  return (*this);
}
const ModelParser& ModelParser::assign(const std::string& tag,
                                       std::string* value,
                                       const std::string& def) const {
  if (this->hasTag(tag)) {
    (*value) = (this->value(tag));
  } else {
    (*value) = def;
  }
  return (*this);
}
