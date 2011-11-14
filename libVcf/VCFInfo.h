#ifndef _VCFINFO_H_
#define _VCFINFO_H_



class VCFInfoValue{
  public:
    int fingerMark; 
    VCFValue* value;
    VCFInfoValue(){
        this->value = new VCFValue;
    }
    ~VCFInfoValue(){
        delete this->value;
        this->value = NULL;
    }
};

class VCFInfo{
public:
    const char* getTag(const char* tag) {
        if (!tag || tag[0] == '\0') 
            return NULL;

        std::string s = tag;
        if (!hasParsed){
            this->parseActual();
        }
        this->tableIter = this->table.find(s);
        if (this->tableIter != this->table.end()){
            return this->tableIter->second->value->toStr();
        } else {
            return NULL;
        }
    } ;
    ~VCFInfo(){
        for (this->tableIter = this->table.begin(); 
             this->tableIter != this->table.end(); 
             this->tableIter ++){
            if (this->tableIter->second != NULL){
                delete this->tableIter->second;
            }
        }
    };
    void reset() { this-> hasParsed = false;};
    void parse(VCFValue* v) {
        this->value = v;
        this->hasParsed = false;
        this->fingerMark ++ ;
    };
    void parseActual(){
        this->hasParsed = true;
        const char* line = this->value->line;
        int b = this->value->beg;
        int e = this->value->beg;
        static std::string key;

        while ( e < this->value->end){
            key.clear();
            // find tag name;
            while(line[e] != '='){
                key.push_back(this->value->line[e++]);
            }
            b = e + 1; // skip '='
        
            this->tableIter = this->table.find(key);
            if ( this->tableIter == this->table.end()){
                VCFInfoValue* f = new VCFInfoValue;
                this->table[key] = f;
                e = parseTillChar(";\t", line, b, f->value);
                f->fingerMark = this->fingerMark;
            } else {
                e = parseTillChar(";\t", line, b, this->tableIter->second->value);
                this->tableIter->second->fingerMark = this->fingerMark;
            };
            e ++ ;
        }
    };
private:
    bool hasParsed;
    VCFValue* value;
    int fingerMark;
    std::map<std::string, VCFInfoValue*> table;
    std::map<std::string, VCFInfoValue*>::iterator tableIter;
}; // VCFInfo

#endif /* _VCFINFO_H_ */
