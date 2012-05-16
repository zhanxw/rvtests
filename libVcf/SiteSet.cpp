#include "SiteSet.h"
#include "IO.h"
#include "TypeConversion.h"

int SiteSet::addSiteFile(const char* fileName){
    int n = 0;
    LineReader lr(fileName);
    std::vector<std::string> fd;
    int pos;
    while(lr.readLineBySep(&fd, " \t")){
        if (fd.size() < 2) continue;
        pos = atoi(fd[1]);
        this->addSite(fd[0], pos);
    };
    return n;
};

// NOTE:
// BED file is [beg, end)
// here we don't care it is 0-based or 1-based and it's up to user to decide
int SiteSet::addBEDFile(const char* fileName){
    int n = 0;
    LineReader lr(fileName);
    std::vector<std::string> fd;
    int beg;
    int end;
    while(lr.readLineBySep(&fd, " \t")){
        if (fd.size() < 3) continue;
        beg = atoi(fd[1]);
        end = atoi(fd[2]);
        for (int i = beg ; i < end; i++) 
            this->addSite(fd[0], i);
    };
    return n;
};
