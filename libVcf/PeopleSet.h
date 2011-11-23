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


#endif /* _PEOPLEINDEX_H_ */
