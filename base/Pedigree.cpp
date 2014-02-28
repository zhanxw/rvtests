#include "Pedigree.h"
#include <assert.h>
#include "Utils.h"

namespace zhanxw{


void Family::addPerson(int pid) {
  people.insert(pid);
  if (this->ped->isFounder(pid))
    founder.insert(pid);
}

void Person::addFather(int id) {
  if (id >= 0) {
    father = (id);
    this->ped->removeFounder(this->id);
    this->founder = false;
  }
};
void Person::addMother(int id) {
  if (id >= 0) {
    mother = (id);
    this->ped->removeFounder(this->id);
    this->founder = false;
  }
}

void Person::addChild(int id){
  if (id >= 0) {
    child.insert(id);
  }
};

/**
 * if @param family or @param person equals to "0", nothing will happen
 * or a family and/or a person will be created if needed
 */
int Pedigree::add(const std::string& family,
                  const std::string& person){
  if (family == "0" || person == "0")
    return 0;

  int fid, pid;
  fid = getFamilyID(family);
  pid = getPersonID(person);

  if (fid < 0 && pid < 0) { // both new
    addFamily(family);
    addPerson(person);
    fid = getFamilyID(family);
    pid = getPersonID(person);
    this->family[fid].addPerson(pid);
    this->people[pid].setFamily(fid);
    return 0;
  } else if (fid < 0 && pid >= 0) { // new fam, old indv,
    fprintf(stderr, "Duplicated person [ %s ] in a new family [ %s ], this is not a supported feature.\n", person.c_str(), family.c_str());
    return -1;
  } else if (fid >= 0 && pid < 0) { // known fam, new indv
    addPerson(person);
    pid = getPersonID(person);
    this->family[fid].addPerson(pid);
    this->people[pid].setFamily(fid);
    return 0;
  } else { // both known
    bool OK = this->family[fid].containPerson(pid)  && this->people[pid].getFamily() == fid;
    if (!OK) {
      fprintf(stderr, "Duplicated person [ %s ] in the duplicated family [ %s ], but not consistent to previous entries\n", person.c_str(), family.c_str());
      return -1;
    }
    return 0;
  }
  return 0;
};

int Pedigree::add(const std::string& family,
                  const std::string& person,
                  const std::string& father,
                  const std::string& mother)  {

  if (person != "0" && (person == father || person == mother) ) {
    fprintf(stderr, "One and his father/mother are both [ %s ]\n", person.c_str());
    return -1;
  }
  if (father != "0" && mother != "0" && father == mother) {
    fprintf(stderr, "Father and mother are both [ %s ]\n", father.c_str());
    return -1;
  }

  // create three pairs
  if (this->add(family, person) < 0 ||
      this->add(family, father) < 0 ||
      this->add(family, mother) < 0) {
    return -1;
  }

  // add relationship
  int fid, pid, fatherId, motherId;
  fid = getFamilyID(family);
  assert(fid >= 0);
  pid = getPersonID(person);
  fatherId = getPersonID(father);
  motherId = getPersonID(mother);

  if (pid >= 0) {
    this->people[pid].addFather(fatherId);
    this->people[pid].addMother(motherId);
  }
  bool errorOccured = false;
  if (fatherId >= 0) {
    this->people[fatherId].addChild(pid);
    if (this->people[fatherId].setGender(MALE) < 0) {
      fprintf(stderr, "Father but not male according to pedigree!\n");
      errorOccured = true;
    };
  }
  if (motherId >= 0) {
    this->people[motherId].addChild(pid);
    if (this->people[motherId].setGender(FEMALE) < 0 ) {
      fprintf(stderr, "Mother but not female according to pedigree!\n");
      errorOccured = true;
    }
  }
  if (errorOccured)
    return -1;
  return 0;
} // Pedigree::add

} // end namespace zhanxw

int loadPedigree(const std::string& fn, zhanxw::Pedigree* ped) {
  zhanxw::Pedigree& p = *ped;
  std::vector<std::string> fd;
  LineReader lr(fn);
  int lineNo = 0;
  while (lr.readLineBySep(&fd, " \t")) {
    ++lineNo;
    removeEmptyField(&fd);
    if (fd.empty()) {
      continue;
    }
    // skip headers
    if (toupper(fd[0]) == "FAM" || toupper(fd[0]) == "FID") {
      continue;
    }
    if (p.add(fd[0], fd[1], fd[2], fd[3]) < 0) {
      fprintf(stderr, "Encounter error when adding line %d.\n", lineNo);
    }
    if (fd.size() >= 5 )
      p.addGender(fd[1], fd[4]);
  }
  
  return 0;
}
