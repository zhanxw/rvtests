#ifndef _VCFINDIVIDUAL_H_
#define _VCFINDIVIDUAL_H_

#include "VCFFunction.h"
#include "VCFValue.h"

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
     *
     * @return 
     */
    void parse() {
        // need to consider missing genotype
        // need to consider missing field
        
        // skip to next 
        if (!this->isInUse()) {
            return;
        }
        
        this->fd.clear();

        VCFValue v;
        int beg = this->self.beg ;

        while(this->self.parseTill(':', beg, &v) == 0) {
            v.line[v.end] = '\0';
            fd.push_back(v);
            beg = v.end + 1;
        };
        /* while (line[v.end] != '\0') { // parse individual values */
        /*     v.end = parseTillChar(":\t\0", line, v.beg, &v); */
        /*     fd.push_back(v); */
        /*     if (line[v.end] == '\t' || line[v.end] == '\0') */
        /*         break; */
        /*     v.beg = v.end + 1; */
        /* } */
        /* this->data.end = v.end; */
        /* return this->data.end; */
        if (fd.size() == 0) {
            fprintf(stderr, "very strange!!");
        }
    };

    const std::string& getName() const {return this->name;};
    void setName(std::string& s) {this->name = s;};
    void include() {this->inUse = true;};
    void exclude() {this->inUse = false;};
    bool isInUse() {return this->inUse;};
    const VCFValue& operator [] (const unsigned int i) const {
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
    VCFValue& getSelf() {return this->self;};
    void rebuildString(std::string* s) {
        assert(s);
        s->clear();
        for (int i = 0; i < fd.size(); i++){
            if (i)
                s->push_back(':');
            *s += this->fd[i].toStr();
        }
    };
  public:
    VCFValue self;            // whole field for the individual
  private:
    bool inUse;
    std::string name;         // id name
    std::vector<VCFValue> fd; // each field separated by ':'
}; // end VCFIndividual

#endif /* _VCFINDIVIDUAL_H_ */
