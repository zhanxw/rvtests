#include "Argument.h"

#define SEP(x) fprintf(stdout, "------------------------- %s -------------------------\n", #x);

int main(int argc, char** argv) {
    BEGIN_PARAMETER_LIST(pl)
        ADD_PARAMETER_GROUP(pl, "Group1")
        ADD_BOOL_PARAMETER(pl, isHelp, "-h", "provide help")
        ADD_INT_PARAMETER(pl, inum, "--number", "a int number")
        ADD_DOUBLE_PARAMETER(pl, dnum, "--double", "a double number")
        ADD_PARAMETER_GROUP(pl, "Group2")
        ADD_STRING_PARAMETER(pl, str, "--string", "sample string")

        ADD_PARAMETER_GROUP(pl, "Input/Output")
        ADD_STRING_PARAMETER(pl,inputVCFFileName, "--input", "input file name")
        ADD_STRING_PARAMETER(pl,outputPrefix, "--outputPrefix", "output prefix for 012 files")

        ADD_PARAMETER_GROUP(pl, "People Filter")
        ADD_STRING_PARAMETER(pl, peopleIncludeID, "--peopleIncludeID", "include people with these ID in the analysis")
        ADD_STRING_PARAMETER(pl, peopleIncludeFile, "--peopleIncludeFile", "specify a file, so ID in that file will be included in analysis")
        ADD_STRING_PARAMETER(pl, peopleExcludeID, "--peopleExcludeID", "exclude people with these ID in the analysis")
        ADD_STRING_PARAMETER(pl, peopleExcludeFile, "--peopleExcludeFile","specify a file, so ID in that file will be excluded in analysis")

        ADD_PARAMETER_GROUP(pl, "INFO field Grepper")
        ADD_STRING_PARAMETER(pl, infoGrep, "--infoGrep", "You can use regular expression to filter the INFO field of the VCF file. "
                                                         "For example, use --infoGrep ANNO=synonymous|ANNO=nonsynonymous to grep both "
                                                         "synonymous and nonsynonymous annotations")
        END_PARAMETER_LIST(pl)
        ;    

    SEP(Begin)
    pl.Read(argc, argv);

    SEP(Status);
    pl.Status();

    SEP(Help);
    pl.Help();

    SEP(Use flags);
    fprintf(stdout, "number = %d\n", FLAG_inum);
    fprintf(stdout, "double = %lf\n", FLAG_dnum);
    fprintf(stdout, "str = %s\n", FLAG_str.c_str());
    fprintf(stdout, "remaining(positional) arguments = ");
    for (int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
        fprintf(stdout, "%s\t", FLAG_REMAIN_ARG[i].c_str());
    }
    fprintf(stdout, "\n");

    SEP(Write to test.param);
    pl.WriteToFile("test.param");
    // you can output test.param by:
    // system("cat test.param");
    
    SEP(Read from test.param);
    pl.ReadFromFile("test.param");
    pl.Status();
    fprintf(stdout, "Remaining arguments: ");
    for (int i = 0; i < FLAG_REMAIN_ARG.size(); i++){
        fprintf(stdout, "%s\t", FLAG_REMAIN_ARG[i].c_str());
    }
    fprintf(stdout, "\n");
    
    SEP(END);
    return 0;
};
