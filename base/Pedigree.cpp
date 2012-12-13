#include "Pedigree.h"
namespace zhanxw{


void Family::addPerson(int pid) {
  people.insert(pid);
  if (this->ped->isFounder(pid))
    founder.insert(pid);
}
void Person::addFather(int id) {
  if (id >= 0) {
    father = (id);
    if (founder) {
      this->ped->removeFounder(id);
    }
    founder = false;
  }
};
void Person::addMother(int id) {
  if (id >= 0) {
    mother = (id);
    if (founder) {
      this->ped->removeFounder(id);
    }
    founder = false;
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
      fprintf(stderr, "Duplicated person [ %s ] in the duplicated family [ %s ], but not consistent to previous entry\n", person.c_str(), family.c_str());
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
  pid = getPersonID(person);
  fatherId = getPersonID(father);
  motherId = getPersonID(mother);

  if (pid >= 0) {
    this->people[pid].addFather(fatherId);
    this->people[pid].addMother(motherId);
  }
  if (fatherId >= 0) {
    this->people[fatherId].addChild(pid);
    this->people[fatherId].setGender(MALE);
  }
  if (motherId >= 0) {
    this->people[motherId].addChild(pid);
    this->people[motherId].setGender(FEMALE);
  }

  return 0;
}


} // end namespace zhanxw

int loadPedigree(const char* fn, zhanxw::Pedigree* ped) {
  zhanxw::Pedigree& p = *ped;
  std::vector<std::string> fd;
  LineReader lr(fn);
  while (lr.readLineBySep(&fd, " \t")) {
    removeEmptyField(&fd);
    if (fd.empty()) {
      continue;
    }
    p.add(fd[0], fd[1], fd[2], fd[3]);
    if (fd.size() >= 5 )
      p.addGender(fd[1], fd[4]);
  };
  return 0;
}
