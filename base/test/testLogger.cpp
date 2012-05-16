#include "Logger.h"

int main(){
    {
        Logger l;
        l.info("hello");
    }
    
    { 
        Logger l("info.log", Logger::INFO);
        l.info("info");
        l.warn("warn");
        l.error("error");
        l.fatal("fatal");
    }
    { 
        Logger l("error.log", Logger::ERROR);
        l.info("info");
        l.warn("warn");
        l.error("error");
        l.fatal("fatal");
    }
};
