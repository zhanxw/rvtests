#ifndef SOCKET_H
#define SOCKET_H

#include <string>

struct addrinfo;

// currently this class is a wrapper to read data in TCP IPv4
class Socket{
 public:
  Socket(const std::string& host,
         int port);
  ~Socket();
  int connect(const std::string& host,
              int port);
  void close();
  int send(const std::string& msg);
  int recv(void* buf, int len);
  bool isUsable() const {return this->usable;}
 private:
  struct addrinfo* servinfo;
  int fd;
  bool usable;
  char buffer[4096];
};

#endif /* SOCKET_H */
