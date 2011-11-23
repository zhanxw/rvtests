#include "RangeList.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>

int main(int argc, char *argv[])
{
    
    RangeList rl;
    rl.addRangeFile("RangeList_test_input");
    std::string range;
    for (unsigned int i = 0; i < rl.size(); i++) {
        rl.obtainRange(i, &range);
        printf("%s\n", range.c_str());
    }
    
    std::string s = "a b\"MID\" c d";
    std::vector<std::string> result;
    unsigned int ret = stringTokenize(s, ' ', &result);
    printf("ret = %d\n", ret);
    dumpStringVector(result);

    ret = stringTokenize(s, "\" ", &result);
    printf("ret = %d\n", ret);
    dumpStringVector(result);

    s = "";
    ret = stringTokenize(s, " ", &result);
    printf("ret = %d\n", ret);
    dumpStringVector(result);

    return 0;
}
