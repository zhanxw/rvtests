#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cassert>

#include "Utils.h"

int main(int argc, char *argv[])
{
    std::string s = "a b\"MID\" c d";
    std::vector<std::string> result;
    unsigned int ret = stringTokenize(s, ' ', &result);
    assert(result.size() == 4);
    assert(result[0] == "a");
    assert(result[1] == "b\"MID\"");
    assert(result[2] == "c");
    assert(result[3] == "d");

    ret = stringTokenize(s, "\" ", &result);
    assert(result.size() == 6);
    assert(result[0] == "a");
    assert(result[1] == "b");
    assert(result[2] == "MID");
    assert(result[3] == "");
    assert(result[4] == "c");
    assert(result[5] == "d");

    s = "";
    ret = stringTokenize(s, " ", &result);
    assert(result.size() == 1);
    assert(result[0] == "");

    return 0;
}
