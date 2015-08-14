#include "Socket.h"

#include <string.h> // memset()
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <unistd.h> // close()

Socket::Socket(const std::string& host, int port):servinfo(NULL), fd(-1), usable(false), quiet(false) {
  connect(host, port);
}

Socket::~Socket() {
  this->close();
}

int Socket::connect(const std::string& host, int port) {
  if (port <= 0)
    return -1;

  // 0. close
  this->close();
  
  // 1. Get IP
  int status;
  struct addrinfo hints;
  this->servinfo = NULL;  // will point to the results
  memset(&hints, 0, sizeof hints);   // make sure the struct is empty
  hints.ai_family = AF_UNSPEC;       // don't care IPv4 or IPv6
  hints.ai_socktype = SOCK_STREAM;   // TCP stream sockets
  hints.ai_flags = AI_PASSIVE;       // fill in my IP for me
  char strPort[128];
  sprintf(strPort, "%d", port);
  if ((status = getaddrinfo(host.c_str(), strPort, &hints, &servinfo)) != 0) {
    if (!this->quiet) {
      fprintf(stderr, "getaddrinfo error: %s\n", gai_strerror(status));
    }
    return -1;
  }

  // 2. get socket
  // use the first servinfo TODO: will try second, third ... if the first fails
  this->fd = ::socket(this->servinfo->ai_family, this->servinfo->ai_socktype, this->servinfo->ai_protocol);
  if (this->fd == -1 ) {
    if (!this->quiet) {
      perror("socket() error");
    }
    return -1;
  }
  
  // 3. connect
  if (::connect(this->fd, this->servinfo->ai_addr, this->servinfo->ai_addrlen) == -1) {
    if (!this->quiet) {
      perror("connect() errror");
    }
    return -1;
  }
  
  usable = true;
  return 0;  
}

void Socket::close() {
  if (this->fd) {
    ::close(fd);
    fd = -1;
  }
  if (servinfo) {
    freeaddrinfo(servinfo);  // free the linked-list
    servinfo = NULL;
  }
  this->usable = false;
}

int Socket::send(const std::string& msg) {
  if (!usable)
    return -1;
  int len = msg.size();
  int nSent = 0;
  while (true) {
    nSent += ::send(this->fd, msg.c_str() + nSent, msg.size(), 0);
    if (nSent == len) {
      break;
    }
    if (nSent < 0) {
      perror("recv() error");
      return -1;
    }
  }
  return nSent;
}

int Socket::recv(void* buf, int len) {
  if (!usable)
    return -1;

  int nRecv = ::recv(this->fd, buf, len, 0);
  return (nRecv);
}
