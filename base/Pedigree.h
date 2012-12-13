#ifndef _PEDIGREE_H_
#define _PEDIGREE_H_

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <set>
#include <string>
#include <vector>
#include "IO.h"
#include "TypeConversion.h"

template<class T1, class T2>
class BiMap {
public:
  T2& get(const T1& t1) {
    return d1[t1];
  }
  T1& get(const T2& t2) {
    return d2[t2];
  }
  T2& operator[] (const T1& t1) {
    return d1[t1];
  }
  T1& operator[] (const T2& t2) {
    return d2[t2];
  }
  void add(const T1& t1, const T2& t2) {
    d1[t1] = t2;
    d2[t2] = t1;
  }
  bool exists(const T1& t1) {
    return d1.count(t1) > 0 ? true : false;
  }
  bool exists(const T2& t2) {
    return d2.count(t2) > 0 ? true : false;
  }
private:
  std::map<T1, T2> d1;
  std::map<T2, T1> d2;
};

class StringIntBiMap {
public:
StringIntBiMap() :defaultString(), defaultInt(-1){};
  int& get(const std::string& t1) {
    return d1[t1];
  }
  std::string& get(const int& t2) {
    return d2[t2];
  }

  const int& operator[] (const std::string& t1) const {
    if (d1.count(t1))
      return d1.find(t1)->second;
    return defaultInt;
  }
  int& operator[] ( const std::string& t1)  {
    return d1[t1];
  }

  const std::string& operator[] (const int& t2) const {
    if (d2.count(t2))
      return d2.find(t2)->second;
    return defaultString;
  }
  std::string& operator[] (const int& t2)  {
    return d2[t2];
  }

  void add(const std::string& t1, const int& t2) {
    d1[t1] = t2;
    d2[t2] = t1;
  }
  bool exists(const std::string& t1) const {
    return d1.count(t1) > 0 ? true : false;
  }
  bool exists(const int& t2) const {
    return d2.count(t2) > 0 ? true : false;
  }
  size_t size() const {
    return d1.size();
  }
private:
  std::string defaultString;
  int defaultInt;
  std::map<std::string, int> d1;
  std::map<int, std::string> d2;
};


namespace zhanxw {
  class Pedigree;

  typedef enum {UNKNOW = 0, MALE = 1, FEMALE = 2} Gender;
  //////////////////////////////////////////////////////////////////////
  // Person
  class Person{
 public:
 Person(Pedigree* p): father(-1), mother(-1), founder(true), gender(UNKNOW), ped(p), family(-1) {};
    void addFather(int id);
    int getFather() const {
      return father;
    }
    void addMother(int id);
    int getMother() const {
      return mother;
    }
    void addChild(int id);
    int setGender(int g) {
      if (g < 0 || g > 2) {
        fprintf(stderr, "Invalid gender\n");
        return -1;
      }
      if (this->gender == UNKNOW)  {
        gender = (Gender)g;
        return 0;
      }
      if ( (int) this->gender != g) {
        fprintf(stderr, "Conflicting gender!");
        return -1;
      }
      return 0;
    };

    int setGender(const std::string& gender) {
      int g = atoi(gender);
      return setGender(g);
    }
    Gender getGender() const {
      return gender;
    }
    void setFamily(int fam) {
      if (fam >= 0) 
        this->family = fam;
    }
    int getFamily() const {
      return this->family;
    }
    void setFounder() {
      this->founder = true;
    }
    bool isFounder() const{
      return founder;
    };
    void dump() const {
      printf("founder = %s, gender = %d, father = %d, mother = %d, ", founder ? "yes" : "no", (int)gender, father, mother);
      printf("child = ");
      for (std::set<int>::const_iterator it = child.begin(); it != child.end(); ++it){
        printf("%d, ", *it);
      }
      puts("");
    };
 private:
    int father;
    int mother;
    std::set<int> child;
    bool founder;
    Gender gender;
    Pedigree* ped;
    int family; // <0: meaning this person is mentioned, but does not have its own entry.
  }; // class Person

  //////////////////////////////////////////////////////////////////////
  // Family
  class Family{
 public:
 Family(Pedigree* p): ped(p) {};
    void addPerson(int pid);
    bool containPerson(int pid) const{
      return (people.count(pid) > 0 ? true : false);
    }
    void removeFounder(int pid) {
      this->founder.erase(pid);
    }
    const std::set<int>& getPeople() const{return people;};
 private:
    std::set<int> people;
    std::set<int> founder;
    Pedigree* ped;
  };

  //////////////////////////////////////////////////////////////////////
  // class Pedigree
  class Pedigree {
 public:
 Pedigree() : totalPeople(0), totalFounder(0){};
  
    int add(const std::string& family,
            const std::string& person);

    int add(const std::string& family,
            const std::string& person,
            const std::string& father,
            const std::string& mother);
    int addGender(const std::string& person, const std::string& gender) {
      int id = getPersonID(person);
      if (id < 0) return -1;
      if (this->people[id].setGender(atoi(gender))) {
        fprintf(stderr, "Failed to set gender [ %s ] for person [ %s ]\n", gender.c_str(), person.c_str());
        return -1;
      }
      return 0;
    };
    /**
     * @return >=0 family id
     */
    int requestFamilyID(const std::string& family) const{
      if (family.empty() || family == "0") return -1;
      int fid;
      // get family id
      if (familyId.exists(family)) {
        fid = getFamilyID(family);
      } else {
        fid = familyId.size();
      }
      return fid;
    }
    /**
     * @return >=0 person id, or -1, if person == "0"
     */
    int requestPersonID(const std::string& person) const{
      if (person.empty() || person == "0") return -1;
      int pid;
      // get person id
      if (personId.exists(person)) {
        pid = getPersonID(person);
      } else {
        pid = personId.size();
      }
      return pid;
    }
    /**
     * For invalid @param id, return true
     * Otherwise, return the true founder status
     */
    bool isFounder(int pid) {
      if (pid < 0 || pid >= (int)family.size())
        return true;
      return people[pid].isFounder();
    }
    int getFamilyID(const std::string& fam) const {
      if (familyId.exists(fam)) {
        return familyId[fam];
      }
      return -1;
    }
    const char* getFamilyName(int fam) const {
      if (familyId.exists(fam)) {
        return familyId[fam].c_str();
      }
      return NULL;
    }

    int getPersonID(const std::string& person) const {
      if (personId.exists(person)) {
        return personId[person];
      }
      return -1;
    }
    const char* getPersonName(int person) const {
      if (personId.exists(person)) {
        return personId[person].c_str();
      }
      return NULL;
    }

    size_t getFamilyNumber() const {
      return family.size();
    }
    const std::vector<Family>& getFamily() const{
      return family;
    }
    size_t getPeopleNumber() const {
      return people.size();
    }
    const std::vector<Person>& getPeople() const {
      return people;
    }
    void addFamily(const std::string& fam) {
      if (familyId.exists(fam)) return;
      Family f(this);
      int n = this->familyId.size();
      this->familyId[fam] = n;
      this->familyId[n] = fam;
      this->family.push_back(f);
    }
    void addPerson(const std::string& person) {
      if (personId.exists(person)) return;
      Person p(this);
      int n = this->personId.size();
      this->personId[person]  = n;
      this->personId[n] = person;
      this->people.push_back(p);
      ++totalPeople;
      ++totalFounder;
    }
    void removeFounder(int id) {
      if (this->people[id].isFounder()) {
        --totalFounder;
        this->family[  this->people[id].getFamily() ] .removeFounder(id);
      }
    }
    int getFounderNumber() const {
      return this->totalFounder;
    }
    /**
     * Put a sequence to @param v such that founders are at the beginning, and all individuals parents are always before itself
    */
    int calculateIterationSequence(std::vector<int>* seq) const {
      std::vector<int>& v = *seq;
      std::set<int> processed;
      v.resize(totalPeople);
      int idx = 0;
      for (size_t i = 0; i < people.size(); ++i) {
        if (people[i].isFounder()) {
          v[idx++] = i;
          processed.insert(i);
        }
      }
      int added = 0;
      while (idx < totalPeople) {
        added = 0;
        for (size_t i = 0; i < people.size(); ++i) {
          if (processed.count(i)) continue;
          if (processed.count ( people[i].getFather() ) &&
              processed.count ( people[i].getMother() )) {
            ++added;
            v[idx++] = i;
            processed.insert(i);
          }
        }
        if (added == 0) {
          fprintf(stderr, "some people are not added, pedigree may have problem\n");
          return -1;
        }
      }
      return 0;
    }
 private:
    std::vector<Family> family;
    std::vector<Person> people;
    StringIntBiMap familyId;
    StringIntBiMap personId;
    int totalPeople;
    int totalFounder;
  }; // class Pedigree

}

int loadPedigree(const char* fn, zhanxw::Pedigree* ped);

inline void dumpPedigree(const zhanxw::Pedigree& ped) {
  int nFam = ped.getFamilyNumber();
  printf("Total %d family loaded: ", nFam);
  for (int i = 0; i < nFam; ++i) {
    printf(" %s, ", ped.getFamilyName(i));
  }
  puts("");

  printf("Total %d founder loaded: \n", ped.getFounderNumber());
  
  int nPeople = ped.getPeopleNumber();
  printf("Total %d people loaded: \n", nPeople);
  for (int i = 0; i < nPeople; ++i) {
    printf("%s: ", ped.getPersonName(i));
    ped.getPeople()[i].dump();
  }
  puts("");

}

inline void dumpPedFile(const zhanxw::Pedigree& ped) {
  const std::vector<zhanxw::Family>& fam = ped.getFamily();
  int nFam = fam.size();
  for (int i = 0; i < nFam; ++i ) {
    const zhanxw::Family & f = fam[i];
    const std::set<int> people = f.getPeople();
    // int nPerson = people.size();
    for (std::set<int>::const_iterator iter = people.begin();
         iter != people.end();
         ++iter) {
      int id = *iter;
      const zhanxw::Person& p = ped.getPeople()[id];
      fprintf(stdout, "%s\t", ped.getFamilyName(i));
      fprintf(stdout, "%s\t", ped.getPersonName(id));
      const char* name = ped.getPersonName(p.getFather());
      fprintf(stdout, "%s\t", name == NULL ? "0" : name);
      name = ped.getPersonName(p.getMother());
      fprintf(stdout, "%s\t", name == NULL ? "0" : name);
      fprintf(stdout, "%d\n", (int)p.getGender());
    }
  }
}


#endif /* _PEDIGREE_H_ */
