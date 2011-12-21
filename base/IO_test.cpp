#include "IO.h"

int main(int argc, char *argv[])
{
    const char fn[] = "README.html.bz2";
    std::string line;
    
    LineReader lr(fn);
    while(lr.readLine(&line)>0) {
        printf("%s\n", line.c_str());
    }

    LineReader lr2(fn);
    std::vector <std::string> fd;
    int lineNo = 0;
    while (lr.readLineBySep(&fd, "\t ")){
        for (int i = 0; i < fd.size(); i++) {
            printf("%s ", fd[i].c_str());
        }
        printf("\n");
    }
    return 0;
}
