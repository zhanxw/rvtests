#include "VCFRecord.h"

//////////////////////////////////////////////////////////////////////
// Code related with include/exclude people
void VCFRecord::includePeople(const std::string& name) {
  if (name.size() == 0) return;

  if (name.find(',') != std::string::npos) {
    std::vector<std::string> v;
    stringTokenize(name, ",", &v);
    includePeople(v);
    return;
  }

  bool included = false;
  for (unsigned int i = 0; i != this->allIndv.size(); i++) {
    VCFIndividual* p = this->allIndv[i];
    if (p->getName() == name) {
      p->include();
      included = true;
    }
  }
  if (!included) {
    fprintf(stderr, "Failed to include sample [ %s ] - not in VCF file.\n",
            name.c_str());
  }
  this->hasAccess = false;
}

void VCFRecord::includePeople(const std::set<std::string>& name) {
  if (name.empty()) return;
  size_t numIncluded = 0;
  for (size_t i = 0; i != this->allIndv.size(); i++) {
    VCFIndividual* p = this->allIndv[i];
    if (name.count(p->getName())) {
      p->include();
      numIncluded++;
    }
  }
  this->hasAccess = false;

  if (numIncluded != name.size()) {
    fprintf(stderr,
            "Warning: Intent to include [ %d ] samples, but actually [ %d ] "
            "samples are included\n",
            (int)name.size(), (int)numIncluded);
  }
}
void VCFRecord::includePeople(const std::vector<std::string>& v) {
  std::set<std::string> s;
  makeSet(v, &s);
  this->includePeople(s);
}
void VCFRecord::includePeopleFromFile(const char* fn) {
  if (!fn || strlen(fn) == 0) return;
  LineReader lr(fn);
  std::vector<std::string> fd;
  std::set<std::string> toInclude;
  while (lr.readLineBySep(&fd, "\t ")) {
    for (unsigned int i = 0; i < fd.size(); i++) toInclude.insert(fd[i]);
  }
  this->includePeople(toInclude);
}
void VCFRecord::includeAllPeople() {
  for (unsigned int i = 0; i != this->allIndv.size(); i++) {
    VCFIndividual* p = this->allIndv[i];
    p->include();
  }
  this->hasAccess = false;
}
void VCFRecord::excludePeople(const std::string& name) {
  if (name.empty()) return;
  for (size_t i = 0; i != this->allIndv.size(); i++) {
    VCFIndividual* p = this->allIndv[i];
    if (p->getName() == name) {
      p->exclude();
    }
  }
  this->hasAccess = false;
}
void VCFRecord::excludePeople(const std::set<std::string>& name) {
  if (name.empty()) return;
  for (size_t i = 0; i != this->allIndv.size(); i++) {
    VCFIndividual* p = this->allIndv[i];
    if (name.count(p->getName())) {
      p->exclude();
    }
  }
  this->hasAccess = false;
}
void VCFRecord::excludePeople(const std::vector<std::string>& v) {
  std::set<std::string> s;
  makeSet(v, &s);
  this->excludePeople(s);
}
void VCFRecord::excludePeopleFromFile(const char* fn) {
  if (!fn || strlen(fn) == 0) return;
  LineReader lr(fn);
  std::vector<std::string> fd;
  std::set<std::string> toExclude;
  while (lr.readLineBySep(&fd, "\t ")) {
    for (unsigned int i = 0; i != fd.size(); i++) toExclude.insert(fd[i]);
  }
  this->excludePeople(toExclude);
}
void VCFRecord::excludeAllPeople() {
  for (unsigned int i = 0; i != this->allIndv.size(); i++) {
    VCFIndividual* p = this->allIndv[i];
    p->exclude();
  }
  this->hasAccess = false;
}
