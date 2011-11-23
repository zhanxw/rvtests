#include "IO.h"

int main(int argc, char *argv[])
{
    const char fn[] = "README.html.bz2";
    std::string line;
    
    LineReader lr(fn);
    while(lr.readLine(&line)>0) {
        printf("%s\n", line.c_str());
    }
    return 0;
}
