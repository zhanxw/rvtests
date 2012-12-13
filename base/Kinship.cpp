#include "Kinship.h"
#include <map>
#include <utility>
namespace zhanxw{
typedef std::map<  std::pair<int, int>, double > Cache;

void updateCache(Cache* cache, int i, int j, double kin) {
  if (i <= j) 
    (*cache)[ std::make_pair(i, j) ] = kin;
  else
    (*cache)[ std::make_pair(j, i) ] = kin;
};

/**
 * @param i th people should be ancestrier than @param j th people
 */
double getKinship(int i, int j, const Pedigree& ped, const std::map<int, int>& order, Cache* cache) {
  // enforce order of i and j
  if (order.find(i)->second > order.find(j)->second ) {
    return getKinship(j, i, ped, order, cache);
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
  if ( i == j) {
    if (ped.getPeople()[i].isFounder()) {
      (*cache)[ std::make_pair(i, j) ] = 0.5;
      return 0.5;
    } else {
      double fatherId = ped.getPeople()[i].getFather();
      double motherId = ped.getPeople()[i].getMother();
      double kin = 0.5 * (1 + getKinship(motherId, fatherId, ped, order, cache));
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
  } else if (ped.getPeople()[i].isFounder() && !ped.getPeople()[j].isFounder()) {
    // i founder; j not founder
    double fatherId = ped.getPeople()[j].getFather();
    double motherId = ped.getPeople()[j].getMother();
    // double k1 = getKinship(i, fatherId, ped, order, cache);
    // double k2 = getKinship(i, motherId, ped, order, cache);
    // printf("k1 = %lf, k2 = %lf\n", k1, k2);
    double kin = 0.5 * (getKinship(i, fatherId, ped, order, cache) + getKinship(i, motherId, ped, order, cache));
    updateCache(cache, i, j, kin);
    return kin;
  } else if (!ped.getPeople()[i].isFounder() && ped.getPeople()[j].isFounder()) {
    // i not founder; j founder
    double fatherId = ped.getPeople()[i].getFather();
    double motherId = ped.getPeople()[i].getMother();
    double kin = 0.5 * (getKinship(fatherId, j, ped, order, cache) + getKinship(motherId, j, ped, order, cache));
    updateCache(cache, i, j, kin);
    return kin;
  } else {
    // i not founder; j not founder
    double fatherId = ped.getPeople()[j].getFather();
    double motherId = ped.getPeople()[j].getMother();
    double kin = 0.5 * (getKinship(i, fatherId, ped, order, cache) + getKinship(i, motherId, ped, order, cache));
    updateCache(cache, i, j, kin);
    return kin;
  }
};
int Kinship::constructFromPedigree(const Pedigree& ped) {
  int nPeople = ped.getPeopleNumber();
  mat.resize(nPeople, nPeople) ;

  Cache cache;

  // founder
  std::vector<int> seq;
  if (ped.calculateIterationSequence(&seq)) {
    return -1;
  }
  // order
  std::map<int, int> order;
  for (size_t i = 0; i < seq.size(); ++i) {
    order[seq[i]] = i;
  }
  
  int nFounder = ped.getFounderNumber();
  for (int i = 0; i < nFounder; ++i ) {
    updateCache(&cache, seq[i], seq[i], 0.5);
    mat[ seq[i] ] [ seq[i] ] = 0.5;
  }

  for (int i = nFounder; i < nPeople; ++i) {
    for (int j = 0; j <= i; ++j ) {
      mat[ seq[i] ][ seq[j] ] = getKinship( seq[i], seq[j] , ped, order, &cache);

      if (j != i) {
        mat[ seq[j] ][ seq[i] ] = mat[ seq[i] ][ seq[j] ];
      }
      // fprintf(stderr, "i = %d, j = %d, kin = %lf\n", seq[i], seq[j], mat[seq[i]][seq[j] ]);
    }
  }
  return 0;
};
int Kinship::constructFromGenotype(const SimpleMatrix& m) {
  return 0;
};
}// end namespacd
