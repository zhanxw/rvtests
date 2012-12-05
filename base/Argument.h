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

#include "OrderedMap.h"
#include "TypeConversion.h"

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
        unsigned int idx = this->flagVec.size();
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
            int n = fscanf(fp, "%s", line);
            if ( n == EOF ) {
                fprintf(stderr, "%s is empty\n.", fileName);
            };
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
    void WriteToFile(const char* fileName) {
        this->WriteToFileWithComment(fileName, "");
    };
    void WriteToFile(FILE* fp) {
        this->WriteToFileWithComment(fp, "");
    };
    
    // write all transated parameters to @param fileName
    // the first line with be @param comment
    // or the information when we write to file.
    void WriteToFileWithComment(const char* fileName, const char* comment) {
        FILE* fp = fopen(fileName, "w");
        this->WriteToFileWithComment(fp, comment);
        fclose(fp);
    }
    void WriteToFileWithComment(FILE* fp, const char* comment) {
        if (strlen(comment) > 0) {
            fprintf(fp, "# %s\n", comment);
        } else {
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
                    getlogin(), hostName, ctime(&tt));
        }
        int numFlagOutputted = 0;
        for (unsigned int i = 0; i != this->flagVec.size(); i++){
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
        for (unsigned int i = 0; i != this->ptrRemainingArg->size(); i++) {
            // separate different flags
            if (numFlagOutputted)
                fprintf(fp, " ");
            // output each remaining argument
            fprintf(fp, " \"%s\"", (*this->ptrRemainingArg)[i].c_str());
        }            
        fprintf(fp, "\n");
    };
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

            bool isLongParam; // I did not make use of this variable.
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
            
            // check if variable flag is a predefined in flagVec or not
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
            this->flagInfoMap[idx].isLongParam = isLongParam;
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
                if (argv[i+1] == '\0') {
                    fprintf(stderr, "WARNING: arg \"%s\" does not have value value\n", argv[i]);
                } else if (argv[i+1][0] == '-' && argv[i+1][1] != '\0') {
                    fprintf(stderr, "WARNING: arg \"%s %s\" does not seem take valid value\n", argv[i], argv[i+1]);
                }
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
  <- group name widt ->   <----------- flag name width ------------------------->
  20          3                    55
*/        
        const int GROUP_WIDTH = 20;
        const int FLAG_WIDTH = 55;
        const char SEP[] = " : ";
        const char EMPTY_SEP[] = "   ";

        for (unsigned int groupIndex = 0; groupIndex < this->groupNameFlagIndexMap.size(); ++groupIndex){
            std::string k;
            std::vector<unsigned int> v;
            this->groupNameFlagIndexMap.at(groupIndex, &k, &v);

            // print group header
            fprintf(stderr, "%*s", GROUP_WIDTH, k.c_str());
            fprintf(stderr, "%s", SEP);

            int availableFlagWidth = FLAG_WIDTH;
            char flagBuffer[FLAG_WIDTH * 4] = "";
            bool firstFlagInLine = true;
            for (unsigned int i = 0; i < v.size(); i++) {
                int idx = v[i];
                std::string flag = this->flagVec[idx];
                FlagInfo &fi = this->flagInfoMap[idx];
                void* data = fi.data;
                int flagWidth = 0;
                
                // print the flag to the flagBuffer
                if (!firstFlagInLine)
                    flagWidth += sprintf(flagBuffer+flagWidth, ", ");
                else 
                    firstFlagInLine = false;
                if (fi.isLongParam) {
                    flagWidth += sprintf(flagBuffer+flagWidth, "--%s [", flag.c_str());
                } else {
                    flagWidth += sprintf(flagBuffer+flagWidth, "-%s [", flag.c_str());
                }
                if (fi.isParsed) {
                    switch(fi.pt){
                    case BOOL_TYPE:
                        if (*(bool*)data)
                            flagWidth += sprintf(flagBuffer+flagWidth, "true");
                        else
                            flagWidth += sprintf(flagBuffer+flagWidth, "false");
                        break;
                    case INT_TYPE:
                        flagWidth += sprintf(flagBuffer+flagWidth, "%d", *(int*)data);
                        break;
                    case DOUBLE_TYPE:
                        flagWidth += sprintf(flagBuffer+flagWidth, "%lf", *(double*)data);
                        break;
                    case STRING_TYPE:
                        flagWidth += sprintf(flagBuffer+flagWidth, "%s", ((std::string*)data)->c_str());
                        break;
                    default:
                        fprintf(stderr, "ERROR: That should be a bug, report to zhanxw@gmail.com");
                        return;
                    }
                }
                flagWidth += sprintf(flagBuffer+flagWidth, "]");

                // if there are not enough spaces and th, we will output a new line
                if (availableFlagWidth < flagWidth && flagWidth < FLAG_WIDTH) {
                    fprintf(stderr, "\n");
                    fprintf(stderr, "%*s", GROUP_WIDTH, " ");
                    fprintf(stderr, "%s", EMPTY_SEP);
                    availableFlagWidth = FLAG_WIDTH;
                    if (flagBuffer[0] == ',') {
                        // print flag
                        fprintf(stderr, "%s", flagBuffer + 2); // 2 is length of ", "
                        availableFlagWidth -= (flagWidth - 2);
                    } else {
                        fprintf(stderr, "%s", flagBuffer); // 2 is length of ", "
                        availableFlagWidth -= flagWidth;
                    }
                    firstFlagInLine = false;
                } else {
                    fprintf(stderr, "%s", flagBuffer); 
                    availableFlagWidth -= flagWidth;
                }
            }
            fprintf(stderr, "\n");
        }
    };
    void Help() {
        const int FLAG_WIDTH = 25;
        const int DOC_WIDTH = 55;
        const char SEP[] = " : ";
        const char EMPTY_SEP[] = "   ";

        for (unsigned int groupIndex = 0; groupIndex < this->groupNameFlagIndexMap.size(); ++groupIndex){
            std::string k;
            std::vector<unsigned int> v;
            this->groupNameFlagIndexMap.at(groupIndex, &k, &v);

            // print group header
            fprintf(stderr, "%s\n", k.c_str());

            for (unsigned int i = 0; i < v.size(); i++) {
                int idx = v[i];
                std::string flag = this->flagVec[idx];
                FlagInfo &fi = this->flagInfoMap[idx];
                // void* data = fi.data;
                // int flagWidth = 0;

                if (fi.isLongParam) {
                    fprintf(stderr, "%*s%s%s", (int)(FLAG_WIDTH - flag.size() ), "--", flag.c_str(), SEP);
                } else {
                    fprintf(stderr, "%*s%s%s", (int)(FLAG_WIDTH - flag.size() ), "-", flag.c_str(), SEP);
                }

                for (unsigned int docIndex = 0; docIndex != fi.doc.size(); docIndex++) {
                    if (docIndex != 0 && docIndex % DOC_WIDTH == 0) {
                        fprintf(stderr, "\n");
                        fprintf(stderr, "%*s%s", FLAG_WIDTH, " ", EMPTY_SEP);
                    }
                    fprintf(stderr, "%c", fi.doc[docIndex]);
                }
                fprintf(stderr, "\n");                    
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
};

#define BEGIN_PARAMETER_LIST(pp)                \
    ParameterParser pp;
#define END_PARAMETER_LIST(pp)                  \
    std::vector<std::string> FLAG_REMAIN_ARG;   \
    pp.AddRemainingArg(&FLAG_REMAIN_ARG);
// pp is an instane of ParameterParser
#define ADD_PARAMETER_GROUP(pp, x)              \
    (pp).AddParameterGroup(x);
#define ADD_BOOL_PARAMETER(pp, x, flag, doc)            \
    bool FLAG_##x;                                      \
    (pp).AddParameter(BOOL_TYPE, &FLAG_##x, flag, doc);                             
#define ADD_INT_PARAMETER(pp, x, flag, doc)             \
    int FLAG_##x;                                       \
    (pp).AddParameter(INT_TYPE, &FLAG_##x, flag, doc);                             
#define ADD_DOUBLE_PARAMETER(pp, x, flag, doc)              \
    double FLAG_##x;                                        \
    (pp).AddParameter(DOUBLE_TYPE, &FLAG_##x, flag, doc);                             
#define ADD_STRING_PARAMETER(pp, x, flag, doc)              \
    std::string FLAG_##x;                                   \
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
