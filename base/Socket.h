#ifndef SOCKET_H
#define SOCKET_H

#include <string>

struct addrinfo;

// currently this class is a wrapper to read data in TCP IPv4
class Socket {
 public:
  Socket(const std::string& host, int port);
  Socket(const std::string& host, int port, double timeoutSeconds);
  virtual ~Socket();
  int connect(const std::string& host, int port);
  int timedConnect(const std::string& host, int port, double seconds);
  void close();
  int send(const std::string& msg);
  int recv(void* buf, int len);
  int timedRecv(void* buf, int len, double seconds);
  bool isUsable() const { return this->usable; }
  void enableQuiet() { this->quiet = true; }
  void disableQuiet() { this->quiet = false; }

 private:
  struct addrinfo* servinfo;
  int fd;
  bool usable;
  // char buffer[4096];
  bool quiet;
};

#endif /* SOCKET_H */
