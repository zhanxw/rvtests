#ifndef _TIME_H_
#define _TIME_H_

#include <string>
#include <ctime>

/**
 * @return a string representing current time, without '\n' at the end
 */
std::string currentTime() {
  time_t t = time(NULL);
  std::string s (ctime(&t));
  s = s.substr(0, s.size() - 1);
  return s;
};



#endif /* _TIME_H_ */
