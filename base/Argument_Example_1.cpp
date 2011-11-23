#include "Argument.h"

int main(int argc, char** argv) {
    BEGIN_PARAMETER_LIST(pl)
        ADD_BOOL_PARAMETER(pl, isHelp, "-h", "provide help")
    END_PARAMETER_LIST(pl)
        ;    
    pl.Read(argc, argv);
    if (FLAG_isHelp)
        fprintf(stdout, "-h is set\n");
    else 
        fprintf(stdout, "-h is NOT set\n");
    return 0;
};
