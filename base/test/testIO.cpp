#include "IO.h"
#include "Utils.h"

#include <stdlib.h>

#include <vector>
#include <string>
#include <cassert>

using std::vector;
using std::string;

char a[]="abc\ndef\nHAH!!DFDSLIJFDSO\nadfa";
char buffer[1024];

int main(int argc, char *argv[])
{
    vector<string> t;
    int ret = stringTokenize(a, '\n', &t);
    assert (ret == 4);

    {
        char fn[] = "abc.txt";
        FileWriter fw(fn);
        fw.write(a);
        fw.close();
        
        LineReader lr(fn);
        std::string ln;
        int i = 0;
        while(lr.readLine(&ln)){
            assert( 0 == strcmp ( t[i].c_str(), ln.c_str()) );
            i ++;
        }
    }

    {
        char fn[] = "abc.txt.gz";
        FileWriter fw(fn);
        fw.write(a);
        fw.close();
        
        LineReader lr(fn);
        std::string ln;
        int i = 0;
        while(lr.readLine(&ln)){
            assert( 0 == strcmp ( t[i].c_str(), ln.c_str()) );
            i ++;
        }
    }

    {
        char fn[] = "abc.txt.bz2";
        FileWriter fw(fn);
        fw.write(a);
        fw.close();
        
        LineReader lr(fn);
        std::string ln;
        int i = 0;
        while(lr.readLine(&ln)){
            assert( 0 == strcmp ( t[i].c_str(), ln.c_str()) );
            i ++;
        }
    }

    {
        char fn[] = "abc.txt.bgzip.gz";
        FileWriter fw(fn, BGZIP);
        fw.write(a);
        fw.close();
        
        LineReader lr(fn);
        std::string ln;
        int i = 0;
        while(lr.readLine(&ln)){
            assert( 0 == strcmp ( t[i].c_str(), ln.c_str()) );
            i ++;
        }

    }

    return 0;
}
