#ifndef _PEOPLESET_H_
#define _PEOPLESET_H_

#include "Utils.h"
#include <string>
#include <set>
#include <vector>

// this class provides which column to filtered (0-based index)
class PeopleSet{
public:
    void readID(const char* allPeopleID);
    void readID(const std::vector<std::string> peopleVec);
    void readIDfromFile(const char* fileName);
    int obtainIDfromFile(const char* fileName, std::vector<std::string>* allID);
    bool contain(const char* name) const {
        std::string s = name;
        return this->contain(s);
    };
    bool contain(const std::string& name) const {
        return this->people.find(name) != this->people.end();
    };
    unsigned int size() const {return this->people.size();};
private:
    std::set <std::string> people;
};

// class PeopleSet
void PeopleSet::readID(const char* allPeopleID){
    if (!strlen(allPeopleID)) return;

    std::vector<std::string> sa;
    stringTokenize(allPeopleID, ",", &sa);
    for (int i = 0; i< sa.size(); i++){
        this->people.insert(sa[i]);
    }
}

void PeopleSet::readID(const std::vector<std::string> peopleVec){
    for (unsigned int i = 0; i < peopleVec.size(); i++ ){
        this->people.insert(peopleVec[i]);        
    }
};

void PeopleSet::readIDfromFile(const char* fileName) {
    if (!strlen(fileName)) return;

    std::vector <std::string> id;
    this->obtainIDfromFile(fileName, &id);
    for (uint32_t i = 0; i < id.size(); i++) {
        this->readID(id[i].c_str());
    }
}

/**
 * File content should be:
 * ID1
 * ID2
 * ID3
 * ...
 * (only the first column will be used)
 */
int  PeopleSet::obtainIDfromFile(const char* fileName, std::vector<std::string>* allID) {
    assert(allID);
    allID->clear();
    std::vector<std::string> fd;
    LineReader lr(fileName);
    while (lr.readLineBySep(&fd, " \t")) {
        (*allID).push_back( fd[0].c_str());
    }
    return 0;
}

#endif /* _PEOPLEINDEX_H_ */
