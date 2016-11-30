#ifndef _ARGUMENT_H_
#define _ARGUMENT_H_

#include <string>
#include <vector>

#include "base/OrderedMap.h"

#define UNUSED(x) ((void)(x))

namespace parameter {

typedef enum PARAMETER_TYPE {
  UNSUPPORTED_TYPE,
  BOOL_TYPE,
  INT_TYPE,
  DOUBLE_TYPE,
  STRING_TYPE
} ParameterType;

struct FlagInfo {
 public:
  FlagInfo()
      : pt(UNSUPPORTED_TYPE), data(NULL), isParsed(false), isLongParam(false){};

 public:
  ParameterType pt;  // the parameter type
  void* data;        // where the parsed value stores
  std::string doc;   // documentation
  bool
      isParsed;  // whether we have parsed the flag and store the value in data.
  bool isLongParam;  // --flag i long parameter, -flag is not
};

/**
 * This class holds all registered flags, parses them and is shared among
 * multiple .cpp files as a singleton
 */
class ParameterParser {
 private:
  // make it a singleton
  explicit ParameterParser();
  // non-copyable
  ParameterParser(const ParameterParser&);
  ParameterParser& operator=(const ParameterParser&);

 public:
  static ParameterParser& getInstance() {
    // static ParameterParser pp;
    // return pp;
    static ParameterParser p;
    return p;
  }

  void AddParameterGroup(const char* name);
  // void AddRemainingArg(void* data) ;
  const std::vector<std::string>& getRemainArg() const {
    return this->FLAG_REMAIN_ARG;
  };
  void AddParameter(ParameterType pt, void* data, const char* flag,
                    const char* doc);
  void InitializeValue(ParameterType pt, void* data, const char* flag);
  // read first uncommented line
  void ReadFromFile(const char* fileName);
  void WriteToFile(FILE* fp);
  void WriteToFile(const char* fileName);
  void WriteToFileWithComment(FILE* fp, const char* comment);
  // write all transated parameters to @param fileName
  // the first line with be @param comment
  // or the information when we write to file.
  void WriteToFileWithComment(const char* fileName, const char* comment);
  /**
   * WriteToStreamWithComment() is EXACTLY the rewrite of
   * WriteToFileWithComment()
   * This repetative work is for output LOG file
   */
  void WriteToStreamWithComment(std::ostream& fout, const char* comment);
  void WriteToStream(std::ostream& fout);
  void Read(int argc, char** argv);
  void Status();
  void Help();

 private:
  // all flags are stored here, "--flag" will store "flag" in flagVec
  std::vector<std::string> flagVec;
  // store flag -> flagInfo
  std::map<unsigned int, FlagInfo> flagInfoMap;
  // store current group name
  std::string currentParameterGroupName;
  // store group and flags (their indices in flagVec)that belongs to that group
  OrderedMap<std::string, std::vector<unsigned int> > groupNameFlagIndexMap;
  // aka. positional argument,
  // those not processed arguments
  std::vector<std::string> FLAG_REMAIN_ARG;
  // std::vector<std::string>* ptrRemainingArg;

  // control how to print help page for parameters
  const static int FLAG_WIDTH = 25;
  const static int DOC_WIDTH = 55;
};  // end class ParameterParser

class ParameterRegister {
 public:
  ParameterRegister(ParameterType pt, void* data, const char* flag,
                    const char* doc) {
    ParameterParser::getInstance().AddParameter(pt, data, flag, doc);
  }
};
class ParameterGroupRegister {
 public:
  ParameterGroupRegister(const char* name) {
    ParameterParser::getInstance().AddParameterGroup(name);
  }
};
}  // end namespace

#define PARAMETER_INSTANCE() (parameter::ParameterParser::getInstance())
#define PARSE_PARAMETER(argc, argv)                 \
  const std::vector<std::string>& FLAG_REMAIN_ARG = \
      PARAMETER_INSTANCE().getRemainArg();          \
  UNUSED(FLAG_REMAIN_ARG);                          \
  PARAMETER_INSTANCE().Read(argc, argv);

// Note: __LINE__ is a macro. It needs to be expanded twice to be concatenated.
// http://stackoverflow.com/questions/1597007/creating-c-macro-with-and-line-token-concatenation-with-positioning-macr
#define TOKENPASTE(x, y) x##y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)
#define ADD_PARAMETER_GROUP(x)                                            \
  namespace parameterGroupRegister {                                      \
  parameter::ParameterGroupRegister TOKENPASTE2(GROUP_REG_, __LINE__)(x); \
  }

#define ADD_PARAMETER(type, typeEnum, value, x, flag, doc)                    \
  namespace parameter {                                                       \
  type FLAG_##x = value;                                                      \
  }                                                                           \
  using parameter::FLAG_##x;                                                  \
  namespace parameterRegister {                                               \
  parameter::ParameterRegister PARAM_REG_##x(typeEnum, &FLAG_##x, flag, doc); \
  }

#define ADD_BOOL_PARAMETER(x, flag, doc) \
  ADD_PARAMETER(bool, parameter::BOOL_TYPE, false, x, flag, doc)
#define ADD_INT_PARAMETER(x, flag, doc) \
  ADD_PARAMETER(int, parameter::INT_TYPE, 0, x, flag, doc)
#define ADD_DOUBLE_PARAMETER(x, flag, doc) \
  ADD_PARAMETER(double, parameter::DOUBLE_TYPE, 0.0, x, flag, doc)
#define ADD_STRING_PARAMETER(x, flag, doc) \
  ADD_PARAMETER(std::string, parameter::STRING_TYPE, "", x, flag, doc)

#define ADD_DEFAULT_INT_PARAMETER(x, def, flag, doc) \
  ADD_PARAMETER(int, parameter::INT_TYPE, def, x, flag, doc)
#define ADD_DEFAULT_DOUBLE_PARAMETER(x, def, flag, doc) \
  ADD_PARAMETER(double, parameter::DOUBLE_TYPE, def, x, flag, doc)
#define ADD_DEFAULT_STRING_PARAMETER(pp, x, def, flag, doc) \
  ADD_PARAMETER(std::string, parameter::STRING_TYPE, def, x, flag, doc)

#define PARAMETER_HELP() PARAMETER_INSTANCE().Help();

#define PARAMETER_STATUS() PARAMETER_INSTANCE().Status();

#define DECLARE_PARAMETER(type, x) \
  namespace parameter {            \
  extern type FLAG_##x;            \
  }                                \
  using parameter::FLAG_##x;
#define DECLARE_BOOL_PARAMETER(x) DECLARE_PARAMETER(bool, x)
#define DECLARE_INT_PARAMETER(x) DECLARE_PARAMETER(int, x)
#define DECLARE_DOUBLE_PARAMETER(x) DECLARE_PARAMETER(double, x)
#define DECLARE_STRING_PARAMETER(x) DECLARE_PARAMETER(std::string, x)

// for compatible reasons, keep them here
#define BEGIN_PARAMETER_LIST()
#define END_PARAMETER_LIST()

#if 0
#define BEGIN_PARAMETER_LIST(pp) ParameterParser* pp;
#define END_PARAMETER_LIST(pp)              \
  std::vector<std::string> FLAG_REMAIN_ARG; \
  pp.AddRemainingArg(&FLAG_REMAIN_ARG);

// pp is an instane of ParameterParser
#define ADD_PARAMETER_GROUP(pp, x) (pp).AddParameterGroup(x);
#define ADD_BOOL_PARAMETER(pp, x, flag, doc) \
  (pp).AddParameter(BOOL_TYPE, &FLAG_##x, flag, doc);

#define ADD_INT_PARAMETER(pp, x, flag, doc) \
  int FLAG_##x = 0;                         \
  (pp).AddParameter(INT_TYPE, &FLAG_##x, flag, doc);
#define ADD_DOUBLE_PARAMETER(pp, x, flag, doc) \
  double FLAG_##x = 0.0;                       \
  (pp).AddParameter(DOUBLE_TYPE, &FLAG_##x, flag, doc);
#define ADD_STRING_PARAMETER(pp, x, flag, doc) \
  std::string FLAG_##x;                        \
  (pp).AddParameter(STRING_TYPE, &FLAG_##x, flag, doc);
// default arguments
#define ADD_DEFAULT_INT_PARAMETER(pp, x, def, flag, doc) \
  int FLAG_##x = def;                                    \
  (pp).AddParameter(INT_TYPE, &FLAG_##x, flag, doc);
#define ADD_DEFAULT_DOUBLE_PARAMETER(pp, x, def, flag, doc) \
  double FLAG_##x = def;                                    \
  (pp).AddParameter(DOUBLE_TYPE, &FLAG_##x, flag, doc);
#define ADD_DEFAULT_STRING_PARAMETER(pp, x, def, flag, doc) \
  std::string FLAG_##x = def;                               \
  (pp).AddParameter(STRING_TYPE, &FLAG_##x, flag, doc);

#endif

////////////////////////////////////////////////////////////////////////////////
// Helper functions
void REQUIRE_STRING_PARAMETER(const std::string& flag, const char* msg);

#endif /* _ARGUMENT_H_ */
