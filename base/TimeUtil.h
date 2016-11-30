#ifndef _TIMEUTIL_H_
#define _TIMEUTIL_H_

#include <ctime>
#include <string>

/**
 * @return a string representing current time, without '\n' at the end
 */
inline std::string currentTime() {
  time_t t = time(NULL);
  std::string s(ctime(&t));
  s = s.substr(0, s.size() - 1);
  return s;
};

#endif /* _TIME_H_ */
