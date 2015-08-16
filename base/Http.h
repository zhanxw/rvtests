#ifndef HTTP_H
#define HTTP_H

#include <string>
#include <vector>

class Socket;
class Http{
 public:
  // Make a connect to @param url
  // @param url must be in the format of:
  // http://domain[:port][/path]
  Http(const std::string& url);
  virtual ~Http();
  int read(std::vector<std::string>* content);
  void enableQuiet();
  void disableQuiet();
 private:
  std::string domain;
  std::string path;
  int port;
  Socket* socket;
  char buf[1024];
  const static int bufSize = 1024;
  std::string proxy;
  int proxyPort;
  std::string request;
  bool quiet;  
};

#endif /* HTTP_H */
