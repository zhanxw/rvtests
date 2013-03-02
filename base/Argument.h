#ifndef _ARGUMENT_H_
#define _ARGUMENT_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> // gethostname
#include <string.h> // strlen
#include <time.h>

#include <string>
#include <map>
#include <queue>
#include <vector>
#include <iostream>

#include "OrderedMap.h"
#include "TypeConversion.h"
#include "LineBreaker.h"

#define UNUSED(x) ((void)(x))

typedef enum PARAMETER_TYPE{
  UNSUPPORTED_TYPE,
  BOOL_TYPE,
  INT_TYPE,
  DOUBLE_TYPE,
  STRING_TYPE
} ParameterType;

struct FlagInfo{
 public:
 FlagInfo():
  pt(UNSUPPORTED_TYPE),
    data(NULL),
    isParsed(false),
    isLongParam(false) {};
 public:
  ParameterType pt;   // the parameter type
  void* data;         // where the parsed value stores
  std::string doc;    // documentation
  bool isParsed;      // whether we have parsed the flag and store the value in data.
  bool isLongParam;   // --flag i long parameter, -flag is not
};

class ParameterParser{
public:
  ParameterParser(){
    this->currentParameterGroupName = "Default Parameter Group";
  }
  void AddParameterGroup(const char* name) {
    // change default group name
    this->currentParameterGroupName = name;
  };
  void AddRemainingArg(void* data) {
    this->ptrRemainingArg = (std::vector<std::string>*) data;
  };
  void AddParameter(ParameterType pt, void* data, const char* flag, const char* doc) {
    std::string f = flag;

    // check parameter validity
    if (!data) {
      fprintf(stderr, "ERROR: Invalid data pointer for flag \"%s\".\n", flag);
      return;
    }
    if (f.size() == 0) {
      fprintf(stderr, "ERROR: Invalid flag name \"%s\"\n", flag);
      return;
    }

    FlagInfo fi;
    // chop flag leading - or --
    if (f.size() > 2 && f[0] == '-' && f[1] == '-') {
      f = f.substr(2);
      fi.isLongParam = true;
    } else if (f.size() > 1 && f[0] == '-') {
      f = f.substr(1);
      fi.isLongParam = false;
    } else {
      fprintf(stderr, "WARNING: Flag \"%s\" is not a valid declaration, try \"-%s\" or \"--%s\".\n", flag, flag, flag);
      return;
    }
    // prevent double existence
    std::vector<std::string>::iterator it;
    for (it = flagVec.begin(); it != flagVec.end(); it++) {
      if (*it == flag){
        fprintf(stderr, "WARNING: Duplicate flag \"%s\" and \"%s\"\n", flag, (*it).c_str());
        return;
      }
    }
    // update internal data
    size_t idx = this->flagVec.size();
    this->flagVec.push_back(f);
    groupNameFlagIndexMap[this->currentParameterGroupName].push_back(idx);
    fi.pt = pt;
    fi.data = data;
    fi.doc = doc;
    flagInfoMap[idx] = fi;

    // initialize data
    switch(pt) {
      case BOOL_TYPE:
        *(bool*)data = false;
        break;
      case INT_TYPE:
        *(int*)data = 0;
        break;
      case DOUBLE_TYPE:
        *(double*)data = 0;
        break;
      case STRING_TYPE:
        *(std::string*)data = "";
        break;
      default:
        fprintf(stderr, "WARNING: Unrecognized parameter type for flag \"%s\"\n", flag);
        return;
    }
  };
  // read first uncommented line
  void ReadFromFile(const char* fileName) {
    FILE* fp = fopen(fileName, "r");
    if (!fp) {
      fprintf(stderr, "ERROR: Cannot open parameter file \"%s\"\n", fileName);
      return ;
    }
    char line[1024];
    int index;
    int len;
    bool isCommentLine = true;
    while (isCommentLine){
      int ret = fscanf(fp, "%s", line);
      UNUSED(ret);
      len = strlen(line);
      if (!len) continue;
      for (index = 0; index < len; index++) {
        if ( line[index] == ' ')
          continue;
        if ( line[index] == '#')
          break;
        isCommentLine = false;
        break;
      }
    }

    int argLen = 10;
    int argc = 1; // we will put each parsed results from argc = 1.
    char** argv = (char**) malloc( sizeof(char*) * argc) ;
    if (!argv) {
      fprintf(stderr, "ERROR: malloc() failed!\n");
      return;
    }
    char prog[] = "DUMMY_PROG";
    argv[0] = prog;
    // parse line[index ... (len-1)] to argc and argv
    std::string arg;
    int numPrefixQuote = 0;
    char quoteChar = '\"';
    for (int i = index; i < len; i++) {
      // if inside the quote
      if (numPrefixQuote > 0) {
        if (line[index] != quoteChar) {
          arg.push_back(line[index]);
          continue;
        } else {
          argv[argc] = (char*) malloc(sizeof(char) * arg.size());
          if (!argv[argc]) {
            fprintf(stderr, "ERROR: malloc() failed!\n");
            return;
          }
          strncpy(argv[argc], arg.c_str(), arg.size());
          arg = "";
        }
      } else {
        // not inside the quote
        if (line[index] == ' ') {
          if (arg.size() == 0)
            continue;
          argv[argc] = (char*) malloc(sizeof(char) * arg.size());
          if (!argv[argc]) {
            fprintf(stderr, "ERROR: malloc() failed!\n");
            return;
          }
          strncpy(argv[argc], arg.c_str(), arg.size());
          arg = "";
        } else {
          arg.push_back(line[index]);
        }
      }
      // allocate new space
      if (argc == argLen) {
        argLen *= 2;
        char** temp = (char**)realloc(argv, sizeof(char*)*argLen);
        if (!temp) {
          fprintf(stderr, "ERROR: realloc() failed!\n");
          return;
        }
        if (argv != temp)
          free(argv);
        argv = temp;
      }
    }
    // read in parameter
    this->Read(argc, argv);

    // clean up memory
    for (int i = 1; i < argc; ++i) {
      if (argv[i])
        free(argv[i]);
    }
    if (argv)
      free(argv);

  };
  void WriteToFile(FILE* fp){
    assert(fp);
    this->WriteToFileWithComment(fp, "");
  };
  void WriteToFile(const char* fileName) {
    this->WriteToFileWithComment(fileName, "");
  };
  void WriteToFileWithComment(FILE* fp, const char* comment) {
    assert(fp);
    if (strlen(comment) > 0) {
      fprintf(fp, "# %s\n", comment);
    } else {
      char username[128] = "unknown_user";
      char* pUsername = getenv("LOGNAME");
      if (pUsername != NULL) {
        strncpy(username, pUsername, 128);
      };
      char hostName[128] = "";
      if (gethostname(hostName, 128) == 0) {
        // succes
      } else {
        // failed
        sprintf(hostName, "Unknown");
      }
      time_t tt = time(0);
      // ctime() will output an extra \n
      fprintf(fp, "# ParameterList created by %s on %s at %s",
              username, hostName, ctime(&tt));
    }
    int numFlagOutputted = 0;
    for (size_t i = 0; i < this->flagVec.size(); i++){
      FlagInfo& fi = this->flagInfoMap[i];
      if (fi.isParsed) {
        // separate different flags
        if (numFlagOutputted)
          fprintf(fp, " ");

        if (fi.isLongParam)
          fprintf(fp, "--");
        else
          fprintf(fp, "-");
        fprintf(fp, "%s", flagVec[i].c_str());
        switch(fi.pt) {
          case BOOL_TYPE:
            break;
          case INT_TYPE:
            fprintf(fp, " %d", *(int*)fi.data);
            break;
          case DOUBLE_TYPE:
            fprintf(fp, " %lf", *(double*)fi.data);
            break;
          case STRING_TYPE:
            fprintf(fp, " \"%s\"", ((std::string*)fi.data)->c_str());
            break;
          default:
            fprintf(stderr, "WARNING: Unrecognized parameter type for flag \"%s\"\n", this->flagVec[i].c_str());
            return;
        }
        numFlagOutputted ++;
      }
    }
    for (size_t i = 0; i < this->ptrRemainingArg->size(); i++) {
      // separate different flags
      if (numFlagOutputted)
        fprintf(fp, " ");
      // output each remaining argument
      fprintf(fp, " \"%s\"", (*this->ptrRemainingArg)[i].c_str());
    }
    fprintf(fp, "\n");
  }
  // write all transated parameters to @param fileName
  // the first line with be @param comment
  // or the information when we write to file.
  void WriteToFileWithComment(const char* fileName, const char* comment) {
    FILE* fp = fopen(fileName, "w");
    assert(fp);
    this->WriteToFileWithComment(fp, comment);
    fclose(fp);
  };
  /**
   * WriteToStreamWithComment() is EXACTLY the rewrite of WriteToFileWithComment()
   * This repetative work is for output LOG file
   */
  void WriteToStreamWithComment(std::ostream& fout, const char* comment) {
    if (strlen(comment) > 0) {
      fout << "# " << comment << "\n";
    } else {
      char username[128] = "unknown_user";
      char* pUsername = getenv("LOGNAME");
      if (pUsername != NULL) {
        strncpy(username, pUsername, 128);
      };
      char hostName[128] = "";
      if (gethostname(hostName, 128) == 0) {
        // succes
      } else {
        // failed
        sprintf(hostName, "Unknown");
      }
      time_t tt = time(0);
      // ctime() will output an extra \n
      fout << "# ParameterList created by " << username << " on "<< hostName << " at " << ctime(&tt);
    }
    int numFlagOutputted = 0;
    for (size_t i = 0; i < this->flagVec.size(); i++){
      FlagInfo& fi = this->flagInfoMap[i];
      if (fi.isParsed) {
        // separate different flags
        if (numFlagOutputted)
          fout << " ";

        if (fi.isLongParam)
          fout << "--";
        else
          fout << "-";
        fout << flagVec[i];
        switch(fi.pt) {
          case BOOL_TYPE:
            break;
          case INT_TYPE:
            fout << " " << *(int*)fi.data;
            break;
          case DOUBLE_TYPE:
            fout << " " << *(double*)fi.data;
            break;
          case STRING_TYPE:
            fout << " \""<< *((std::string*)fi.data) <<  "\"";
            break;
          default:
            fprintf(stderr, "WARNING: Unrecognized parameter type for flag \"%s\"\n", this->flagVec[i].c_str());
            return;
        }
        numFlagOutputted ++;
      }
    }
    for (size_t i = 0; i < this->ptrRemainingArg->size(); i++) {
      // separate different flags
      if (numFlagOutputted)
        fout << " ";
      // output each remaining argument
      fout << " \"" << (*this->ptrRemainingArg)[i] << "\"";
    }
    fout << "\n";
  }
  void WriteToStream(std::ostream& fout) {
    this->WriteToStreamWithComment(fout, "");
  }
  void Read(int argc, char** argv) {
    std::string flag;
    for (int i = 1; i < argc; i++) {
      // removing leading - or -- for flags
      flag = argv[i];
      // flag = '--' meaning the end of all parameter
      if (flag == "--") {
        ++i ;
        while (i < argc)
          this->ptrRemainingArg->push_back(argv[i++]);
        return;
      }

      // flag = '-' meaning standard input, we will put it in the remaingArg
      if (flag == "-") {
        this->ptrRemainingArg->push_back(argv[i]);
        continue;
      }

      bool isLongParam = false; // I did not make use of this variable.
      UNUSED(isLongParam);
      unsigned int choppedLeadingDash = 0;
      // user may input ---flag, ----flag, ..., and we will chop all leading -
      // user may also input ---, ----, but I don't understand what that means, so report error
      while (flag.size() > 0){
        if (flag[0] == '-') {
          choppedLeadingDash ++;
          flag = flag.substr(1);
        } else {
          break;
        }
      }
      if (flag.size() > 0) {
        if (choppedLeadingDash > 1) {
          isLongParam = true;
        } else {
          isLongParam = false;
        }
      } else { // (flag.size() == 0
        fprintf(stderr, "ERROR: we don't understand the argument \"%s\"\n", argv[i]);
        return;
      }

      // check if variable flag is a predefined flag or not
      std::vector<std::string>::iterator it;
      for (it = flagVec.begin(); it != flagVec.end(); it++) {
        if (*it == flag){
          break;
        }
      }
      // variable flag will be added to remainingArg
      if (it == this->flagVec.end()) {
        this->ptrRemainingArg->push_back(argv[i]);
        continue;
      }
      // NOTE:
      // we could check if user misuse -- and -
      // but i did not see much use here.

      // parse data
      unsigned int idx = it - this->flagVec.begin();
      void* data = this->flagInfoMap[idx].data;
      if (this->flagInfoMap[idx].isParsed) {
        fprintf(stderr, "WARNING: flag \"%s\" provided more than once, the previous value will be overwritten\n", argv[i]);
      }
      this->flagInfoMap[idx].isParsed = true;
      switch(this->flagInfoMap[idx].pt) {
        case BOOL_TYPE:
          *(bool*)data = true;
          break;
        case INT_TYPE:
          if (!str2int(argv[i+1], (int*)data)) {
            fprintf(stderr, "WARNING: arg \"%s %s\" does not give valid integer\n", argv[i], argv[i+1]);
          }
          ++i;
          break;
        case DOUBLE_TYPE:
          if (!str2double(argv[i+1], (double*)data)) {
            fprintf(stderr, "WARNING: arg \"%s %s\" does not give valid double\n", argv[i], argv[i+1]);
          }
          ++i;
          break;
        case STRING_TYPE:
          *(std::string*)data = argv[++i];
          break;
        default:
          fprintf(stderr, "ERROR: Unrecognized parameter type for flag %s\n", argv[i]);
          return;
      }
    }
  };
  void Status() {
    fprintf(stderr, "The following parameters are available.  Ones with \"[]\" are in effect:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Available Options\n");
    /*
      Format illustration:
      Individual Filter : --indvDepthMin [], --indvDepthMax [], --indvQualMin []
      <- group name widt -><- SEP -><----------- flag name width ------------------------->
      GROUP_WIDTH                               FLAG_WIDTH
    */

    LineBreaker leftColumn(FLAG_WIDTH);
    LineBreaker rightColumn(DOC_WIDTH);
    rightColumn.setSeparator(",");
    std::string left;
    std::string right;
    std::string rightItem;

    for (size_t groupIndex = 0; groupIndex < this->groupNameFlagIndexMap.size(); ++groupIndex){
      std::string k;
      std::vector<unsigned int> v;
      this->groupNameFlagIndexMap.at(groupIndex, &k, &v);

      // process group header
      left = k;
      left += " :";
      right.clear();

      for (size_t i = 0; i < v.size(); i++) { // loop flags under the same param group
        int idx = v[i];
        const std::string& flag = this->flagVec[idx];
        FlagInfo &fi = this->flagInfoMap[idx];
        void* data = fi.data;

        if (fi.isLongParam) {
          rightItem += " --";
          rightItem += flag;
          if (fi.pt != BOOL_TYPE || fi.isParsed) {
            rightItem += " [";
          }
        } else {
          rightItem += " -";
          rightItem += flag;
          if (fi.pt != BOOL_TYPE || fi.isParsed) {
            rightItem += " [";
          }
        }
        if (fi.isParsed) {
          switch(fi.pt){
            case BOOL_TYPE:
              if (*(bool*)data) {
                rightItem += "true";
              } else {
                rightItem += "false";
              }
              break;
            case INT_TYPE:
              rightItem += toString(*(int*)data);
              break;
            case DOUBLE_TYPE:
              rightItem += toString(*(double*)data);
              break;
            case STRING_TYPE:
              rightItem += *(std::string*)data;
              break;
            default:
              fprintf(stderr, "ERROR: That should be a bug, report to zhanxw@gmail.com");
              return;
          }
        }
        if (fi.pt != BOOL_TYPE || fi.isParsed) {
          rightItem += "]";
        }

        right += rightItem;
        right += ",";
        rightItem.clear();
      } // end loop per param group
      // fprintf(stderr, "\n");
      leftColumn.setContent(left);
      rightColumn.setContent(right);
      printTwoColumn(stderr, leftColumn, rightColumn, "");
    }
  };
  void Help() {
    LineBreaker leftColumn(FLAG_WIDTH);
    LineBreaker rightColumn(DOC_WIDTH);
    std::string left;
    for (size_t groupIndex = 0; groupIndex < this->groupNameFlagIndexMap.size(); ++groupIndex){
      std::string k;
      std::vector<unsigned int> v;
      this->groupNameFlagIndexMap.at(groupIndex, &k, &v);

      // print group header
      fprintf(stderr, "%s\n", k.c_str());

      for (size_t i = 0; i < v.size(); i++) {
        int idx = v[i];
        const std::string& flag = this->flagVec[idx];
        FlagInfo &fi = this->flagInfoMap[idx];

        if (fi.isLongParam) {
          left = "--";
          left += flag;
        } else {
          left = "-";
          left += flag;
        }
        left += " :";

        leftColumn.setContent(left);
        rightColumn.setContent(fi.doc);
        printTwoColumn(stderr, leftColumn, rightColumn, " ");
      }
    }
  };
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
  // those not processed argument
  std::vector<std::string>* ptrRemainingArg;
  const static int FLAG_WIDTH = 25;
  const static int DOC_WIDTH = 55;
};

#define BEGIN_PARAMETER_LIST(pp)                \
  ParameterParser pp;
#define END_PARAMETER_LIST(pp)                  \
  std::vector<std::string> FLAG_REMAIN_ARG;     \
  pp.AddRemainingArg(&FLAG_REMAIN_ARG);
// pp is an instane of ParameterParser
#define ADD_PARAMETER_GROUP(pp, x)              \
  (pp).AddParameterGroup(x);
#define ADD_BOOL_PARAMETER(pp, x, flag, doc)            \
  bool FLAG_##x;                                        \
  (pp).AddParameter(BOOL_TYPE, &FLAG_##x, flag, doc);
#define ADD_INT_PARAMETER(pp, x, flag, doc)             \
  int FLAG_##x;                                         \
  (pp).AddParameter(INT_TYPE, &FLAG_##x, flag, doc);
#define ADD_DOUBLE_PARAMETER(pp, x, flag, doc)          \
  double FLAG_##x;                                      \
  (pp).AddParameter(DOUBLE_TYPE, &FLAG_##x, flag, doc);
#define ADD_STRING_PARAMETER(pp, x, flag, doc)          \
  std::string FLAG_##x;                                 \
  (pp).AddParameter(STRING_TYPE, &FLAG_##x, flag, doc);


////////////////////////////////////////////////////////////////////////////////
// Helper functions
static void REQUIRE_STRING_PARAMETER(const std::string& flag, const char* msg){
  if (flag.size() == 0){
    fprintf(stderr, "%s\n", msg);
    abort();
  }
};

#endif /* _ARGUMENT_H_ */
