#include "Kinship.h"
#include <map>
#include <utility>
namespace zhanxw {
typedef std::map<std::pair<int, int>, double> Cache;
void printCache(Cache* c) {
  std::map<std::pair<int, int>, double>::const_iterator it;
  for (it = c->begin(); it != c->end(); ++it) {
    fprintf(stderr, "(%d, %d) = %g\n", it->first.first, it->first.second,
            it->second);
  }
}

void updateCache(Cache* cache, int i, int j, double kin) {
  if (i <= j)
    (*cache)[std::make_pair(i, j)] = kin;
  else
    (*cache)[std::make_pair(j, i)] = kin;
};

/**
 * @param i th people should be ancestrier than @param j th people
 */
double getKinshipFromOrderedPair(int i, int j, const Pedigree& ped,
                                 const std::map<int, int>& order,
                                 Cache* cache) {
  // enforce order of i and j so that order[i] <= order[j]
  if (order.find(i)->second > order.find(j)->second) {
    return getKinshipFromOrderedPair(j, i, ped, order, cache);
  }

  // different family
  if (ped.getPeople()[i].getFamily() != ped.getPeople()[j].getFamily())
    return 0.0;

  // check cache
  Cache::const_iterator iter = cache->find(std::make_pair(i, j));
  if (iter != cache->end()) {
    return iter->second;
  }

  // same family
  // calculate kinship for the same one.
  if (i == j) {
    if (ped.getPeople()[i].isFounder()) {
      (*cache)[std::make_pair(i, j)] = 0.5;
      return 0.5;
    } else {
      double fatherId = ped.getPeople()[i].getFather();
      double motherId = ped.getPeople()[i].getMother();
      double kin = 0.5 * (1. + getKinshipFromOrderedPair(motherId, fatherId,
                                                         ped, order, cache));
      // (*cache)[ std::make_pair(i, j) ] = kin;
      updateCache(cache, i, j, kin);
      return kin;
    }
  }

  //  i != j
  if (ped.getPeople()[i].isFounder() && ped.getPeople()[j].isFounder()) {
    // i founder; j founder
    updateCache(cache, i, j, 0.0);
    return 0.0;
  } else if (ped.getPeople()[i].isFounder() &&
             !ped.getPeople()[j].isFounder()) {
    // i founder; j not founder
    double fatherId = ped.getPeople()[j].getFather();
    double motherId = ped.getPeople()[j].getMother();
    // double k1 = getKinship(i, fatherId, ped, order, cache);
    // double k2 = getKinship(i, motherId, ped, order, cache);
    // printf("k1 = %lf, k2 = %lf\n", k1, k2);
    double kin =
        0.5 * (getKinshipFromOrderedPair(i, fatherId, ped, order, cache) +
               getKinshipFromOrderedPair(i, motherId, ped, order, cache));
    updateCache(cache, i, j, kin);
    return kin;
  } else if (!ped.getPeople()[i].isFounder() &&
             ped.getPeople()[j].isFounder()) {
    // i not founder; j founder
    double fatherId = ped.getPeople()[i].getFather();
    double motherId = ped.getPeople()[i].getMother();
    double kin =
        0.5 * (getKinshipFromOrderedPair(fatherId, j, ped, order, cache) +
               getKinshipFromOrderedPair(motherId, j, ped, order, cache));
    updateCache(cache, i, j, kin);
    return kin;
  } else {
    // i not founder; j not founder
    double fatherId = ped.getPeople()[j].getFather();
    double motherId = ped.getPeople()[j].getMother();
    double kin =
        0.5 * (getKinshipFromOrderedPair(i, fatherId, ped, order, cache) +
               getKinshipFromOrderedPair(i, motherId, ped, order, cache));
    updateCache(cache, i, j, kin);
    return kin;
  }
};

int Kinship::constructFromPedigree(const Pedigree& ped) {
  int nPeople = ped.getPeopleNumber();
  mat.resize(nPeople, nPeople);

  Cache cache;

  // founder
  std::vector<int> seq;
  if (ped.calculateIterationSequence(&seq)) {
    return -1;
  }

  // order
  // e.g. seq = { 2, 4, 1, 3, 5}
  // then order = {2:0, 4:1, 1:2, 3:3, 5:4}
  // higher order meaning less likely to be a founder
  std::map<int, int> order;
  for (size_t i = 0; i < seq.size(); ++i) {
    order[seq[i]] = i;
  }

  int nFounder = ped.getFounderNumber();
  for (int i = 0; i < nFounder; ++i) {
    updateCache(&cache, seq[i], seq[i], 0.5);
    mat[seq[i]][seq[i]] = 0.5;
  }

  for (int i = nFounder; i < nPeople; ++i) {
    for (int j = 0; j <= i; ++j) {
      mat[seq[i]][seq[j]] =
          getKinshipFromOrderedPair(seq[i], seq[j], ped, order, &cache);

      if (j != i) {
        mat[seq[j]][seq[i]] = mat[seq[i]][seq[j]];
      }
      // fprintf(stderr, "i = %d, j = %d, kin = %lf\n", seq[i], seq[j],
      // mat[seq[i]][seq[j] ]);
    }
  }
  return 0;
};

//////////////////////////////////////////////////
// Kinship calculation for sex chromosome
//////////////////////////////////////////////////
/**
 * @param i th people should be ancestrier than @param j th people
 */
double getKinshipFromOrderedPairForX(int i, int j, const Pedigree& ped,
                                     const std::map<int, int>& order,
                                     Cache* cache) {
  // enforce order of i and j so that order[i] <= order[j]
  if (order.find(i)->second > order.find(j)->second) {
    return getKinshipFromOrderedPairForX(j, i, ped, order, cache);
  }

  // different family
  if (ped.getPeople()[i].getFamily() != ped.getPeople()[j].getFamily())
    return 0.0;

  // check cache
  Cache::const_iterator iter = cache->find(std::make_pair(i, j));
  if (iter != cache->end()) {
    return iter->second;
  }

  // same family
  // calculate kinship for the same one.
  if (i == j) {
    // i/j is male
    if (ped.getPeople()[i].getGender() == MALE) {
      (*cache)[std::make_pair(i, j)] = 1.0;
      return 1.0;
    }

    // i/j is female
    if (ped.getPeople()[i].isFounder()) {  // founder female.
      (*cache)[std::make_pair(i, j)] = 0.5;
      return 0.5;
    } else {  // non-founder female
      double fatherId = ped.getPeople()[i].getFather();
      double motherId = ped.getPeople()[i].getMother();
      double kin = 0.5 * (1. + getKinshipFromOrderedPairForX(
                                   motherId, fatherId, ped, order, cache));
      // (*cache)[ std::make_pair(i, j) ] = kin;
      updateCache(cache, i, j, kin);
      return kin;
    }
  }

  //  i != j
  if (ped.getPeople()[i].isFounder() && ped.getPeople()[j].isFounder()) {
    // i founder; j founder
    updateCache(cache, i, j, 0.0);
    return 0.0;
  } else if (ped.getPeople()[i].isFounder() &&
             !ped.getPeople()[j].isFounder()) {
    // i founder; j not founder
    double fatherId = ped.getPeople()[j].getFather();
    double motherId = ped.getPeople()[j].getMother();
    double kin;
    if (ped.getPeople()[j].getGender() == MALE) {
      kin = getKinshipFromOrderedPairForX(i, motherId, ped, order, cache);
    } else {
      kin =
          0.5 * (getKinshipFromOrderedPairForX(i, fatherId, ped, order, cache) +
                 getKinshipFromOrderedPairForX(i, motherId, ped, order, cache));
    }
    updateCache(cache, i, j, kin);
    return kin;
  } else if (!ped.getPeople()[i].isFounder() &&
             ped.getPeople()[j].isFounder()) {
    // i not founder; j founder
    double fatherId = ped.getPeople()[i].getFather();
    double motherId = ped.getPeople()[i].getMother();
    double kin;
    if (ped.getPeople()[i].getGender() == MALE) {
      kin = getKinshipFromOrderedPairForX(motherId, j, ped, order, cache);
    } else {
      kin =
          0.5 * (getKinshipFromOrderedPairForX(fatherId, j, ped, order, cache) +
                 getKinshipFromOrderedPairForX(motherId, j, ped, order, cache));
    }
    updateCache(cache, i, j, kin);
    return kin;
  } else {
    // i not founder; j not founder (and i is not a descendant of j)
    double fatherId = ped.getPeople()[j].getFather();
    double motherId = ped.getPeople()[j].getMother();
    double kin;
    if (ped.getPeople()[j].getGender() == MALE) {
      kin = getKinshipFromOrderedPairForX(i, motherId, ped, order, cache);
    } else {
      kin =
          0.5 * (getKinshipFromOrderedPairForX(i, fatherId, ped, order, cache) +
                 getKinshipFromOrderedPairForX(i, motherId, ped, order, cache));
    }
    updateCache(cache, i, j, kin);
    return kin;
  }
};

int KinshipForX::constructFromPedigree(const Pedigree& ped) {
  int nPeople = ped.getPeopleNumber();
  // make sure sex is known for everyone.
  // or we cannot calculate kinship (e.g. 1 for outbred male, 0.5 for oubred
  // female)
  for (int i = 0; i < nPeople; ++i) {
    if (ped.getPeople()[i].getGender() != MALE &&
        ped.getPeople()[i].getGender() != FEMALE) {
      return -1;
    }
  }

  mat.resize(nPeople, nPeople);

  Cache cache;

  // founder
  std::vector<int> seq;
  if (ped.calculateIterationSequence(&seq)) {
    return -1;
  }

  // order
  // e.g. seq = { 2, 4, 1, 3, 5}
  // then order = {2:0, 4:1, 1:2, 3:3, 5:4}
  // higher order meaning less likely to be a founder
  std::map<int, int> order;
  for (size_t i = 0; i < seq.size(); ++i) {
    order[seq[i]] = i;
  }

  int nFounder = ped.getFounderNumber();
  for (int i = 0; i < nFounder; ++i) {
    if (ped.getPeople()[seq[i]].getGender() == MALE) {
      updateCache(&cache, seq[i], seq[i], 1.0);
      mat[seq[i]][seq[i]] = 1.0;
    } else {
      updateCache(&cache, seq[i], seq[i], 0.5);
      mat[seq[i]][seq[i]] = 0.5;
    }
  }

  for (int i = nFounder; i < nPeople; ++i) {
    for (int j = 0; j <= i; ++j) {
      mat[seq[i]][seq[j]] =
          getKinshipFromOrderedPairForX(seq[i], seq[j], ped, order, &cache);

      if (j != i) {
        mat[seq[j]][seq[i]] = mat[seq[i]][seq[j]];
      }
      // fprintf(stderr, "i = %d, j = %d, kin = %lf\n", seq[i], seq[j],
      // mat[seq[i]][seq[j] ]);
    }
  }
  return 0;
};

}  // end namespacd
