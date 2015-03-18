#include "Regex.h"
#include <cassert>

int main(int argc, char *argv[])
{
    Regex re;
    assert (re.readPattern("") == 0);
    assert(re.match("") == false);

    re.readPattern("abc");
    assert(re.match("abc") == true);

    re.readPattern("abc");
    assert(re.match("__abc") == true);

    re.readPattern("abc");
    assert(re.match("abc__") == true);

    re.readPattern("abc");
    assert(re.match("__abc__") == true);

    re.readPattern("abc");
    assert(re.match("abc___", 2, 4) == false);
    
    re.readPattern("abc|def");
    assert(re.match("abc") == true);

    re.readPattern("abc|def");
    assert(re.match("__def__") == true);

    re.readPattern("abc|def");
    assert(re.match("abef") == false);

    return 0;
}
