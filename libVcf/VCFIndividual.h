#ifndef _VCFINDIVIDUAL_H_
#define _VCFINDIVIDUAL_H_

#include "VCFFunction.h"
#include "VCFValue.h"

#if 0
extern int parseTillChar(const char c, const char* line, const int beg, VCFValue* ret); 
extern int parseTillChar(const char* c, const char* line, const int beg, VCFValue* ret);
#endif

// we assume format are always  GT:DP:GQ:GL
class VCFIndividual{
public:

    // FUNC parseFunction[4];
    VCFIndividual(){
        this->include();  // by default, enable everyone
    };
    /**
     * 0-base index for beg and end, e.g.
     *     0 1 2  3
     *     A B C \t
     * beg = 0, end = 3 (line[end] = '\t' or line[end] = '\0')
     */
    int parse(const char* line, const int beg) {
        // need to consider missing genotype
        // need to consider missing field
        
        // skip to next 
        if (!this->isInUse()) {
            return parseTillChar('\t', line, beg, & this->data);
        }
        
        this->data.line = line;
        this->data.beg = beg;
        
        this->fd.clear();
        VCFValue v;
        v.line = line;
        v.beg = beg;
        v.end = beg;
        while (line[v.end] != '\0') { // parse individual values
            v.end = parseTillChar(":\t\0", line, v.beg, &v);
            fd.push_back(v);
            if (line[v.end] == '\t' || line[v.end] == '\0')
                break;
            v.beg = v.end + 1;
        }
        this->data.end = v.end;
        return this->data.end;
    };

    const std::string& getName() const {return this->name;};
    void setName(std::string& s) {this->name = s;};
    void include() {this->inUse = true;};
    void exclude() {this->inUse = false;};
    bool isInUse() {return this->inUse;};
    const VCFValue& operator [] (const unsigned int i) const __attribute__ ((deprecated)) {
        if (i >= fd.size()){
            FATAL("index out of bound!");
        }
        return (this->fd[i]);                   
    };
    VCFValue& operator [] (const unsigned int i) {
        if (i >= fd.size()){
            FATAL("index out of bound!");
        }
        return (this->fd[i]);
    };
    VCFValue& getData() {return this->data;};
private:
    bool inUse;
    std::string name;

    VCFValue data;            // whole field for the individual
    std::vector<VCFValue> fd; // each field separated by ':'
}; // end VCFIndividual

#endif /* _VCFINDIVIDUAL_H_ */
