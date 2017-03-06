#include "Argument.h"

#define SEP(x)                                                                \
  fprintf(stdout, "------------------------- %s -------------------------\n", \
          #x);

// BEGIN_PARAMETER_LIST(pl)
ADD_PARAMETER_GROUP("Group1")
ADD_BOOL_PARAMETER(isHelp, "-h", "provide help")
ADD_INT_PARAMETER(inum, "--number", "a int number")
ADD_DOUBLE_PARAMETER(dnum, "--double", "a double number")
ADD_PARAMETER_GROUP("Group2")
ADD_STRING_PARAMETER(str, "--string", "sample string")

ADD_PARAMETER_GROUP("Input/Output")
ADD_STRING_PARAMETER(inputVCFFileName, "--input", "input file name")
ADD_STRING_PARAMETER(outputPrefix, "--outputPrefix",
                     "output prefix for 012 files")

ADD_PARAMETER_GROUP("People Filter")
ADD_STRING_PARAMETER(peopleIncludeID, "--peopleIncludeID",
                     "include people with these ID in the analysis")
ADD_STRING_PARAMETER(
    peopleIncludeFile, "--peopleIncludeFile",
    "specify a file, so ID in that file will be included in analysis")
ADD_STRING_PARAMETER(peopleExcludeID, "--peopleExcludeID",
                     "exclude people with these ID in the analysis")
ADD_STRING_PARAMETER(
    peopleExcludeFile, "--peopleExcludeFile",
    "specify a file, so ID in that file will be excluded in analysis")

ADD_PARAMETER_GROUP("INFO field Grepper")
ADD_STRING_PARAMETER(infoGrep, "--infoGrep",
                     "You can use regular expression to filter the INFO "
                     "field of the VCF file. "
                     "For example, use --infoGrep "
                     "ANNO=synonymous|ANNO=nonsynonymous to grep both "
                     "synonymous and nonsynonymous annotations")
// END_PARAMETER_LIST(pl);

int main(int argc, char** argv) {
  SEP(Begin);
  // pl.Read(argc, argv);
  PARSE_PARAMETER(argc, argv);

  SEP(Status);
  // pl.Status();
  PARAMETER_STATUS();

  SEP(Help);
  // pl.Help();
  PARAMETER_HELP();

  SEP(Use flags);
  fprintf(stdout, "number = %d\n", FLAG_inum);
  fprintf(stdout, "double = %lf\n", FLAG_dnum);
  fprintf(stdout, "str = %s\n", FLAG_str.c_str());
  fprintf(stdout, "remaining(positional) arguments = ");
  for (size_t i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
    fprintf(stdout, "%s\t", FLAG_REMAIN_ARG[i].c_str());
  }
  fprintf(stdout, "\n");

  SEP(Write to test.param);
  PARAMETER_INSTANCE().WriteToFile("test.param");
  // you can output test.param by:
  // system("cat test.param");

  SEP(Read from test.param);
  PARAMETER_INSTANCE().ReadFromFile("test.param");
  PARAMETER_INSTANCE().Status();
  fprintf(stdout, "Remaining arguments: ");
  for (size_t i = 0; i < FLAG_REMAIN_ARG.size(); i++) {
    fprintf(stdout, "%s\t", FLAG_REMAIN_ARG[i].c_str());
  }
  fprintf(stdout, "\n");

  SEP(END);
  return 0;
};
