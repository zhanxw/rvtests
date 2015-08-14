#ifndef HTTP_H
#define HTTP_H

#include <string>
#include <vector>

class Socket;
class Http{
 public:
  Http(const std::string& url, int port = 80);
  virtual ~Http();
  int read(std::vector<std::string>* content);
 private:
  std::string domain;
  std::string path;
  Socket* socket;
  char buf[1024];
  const static int bufSize = 1024;
};

#endif /* HTTP_H */
