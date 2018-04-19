#include "Socket.h"

#include <fcntl.h>
#include <netdb.h>
#include <string.h>  // memset()
#include <sys/socket.h>
#include <sys/time.h>  // struct timeval
#include <sys/types.h>
#include <unistd.h>  // close()

Socket::Socket(const std::string& host, int port)
    : servinfo(NULL), fd(-1), usable(false), quiet(false) {
  connect(host, port);
}

Socket::Socket(const std::string& host, int port, double timeoutSeconds)
    : servinfo(NULL), fd(-1), usable(false), quiet(false) {
  timedConnect(host, port, timeoutSeconds);
}

Socket::~Socket() { this->close(); }

int Socket::connect(const std::string& host, int port) {
  if (host.empty()) return -1;
  if (port <= 0) return -1;

  // 0. close
  this->close();

  // 1. Get IP
  int status;
  struct addrinfo hints;
  this->servinfo = NULL;            // will point to the results
  memset(&hints, 0, sizeof hints);  // make sure the struct is empty
  hints.ai_family = AF_UNSPEC;      // don't care IPv4 or IPv6
  hints.ai_socktype = SOCK_STREAM;  // TCP stream sockets
  hints.ai_flags = AI_PASSIVE;      // fill in my IP for me
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
  this->fd = ::socket(this->servinfo->ai_family, this->servinfo->ai_socktype,
                      this->servinfo->ai_protocol);
  if (this->fd == -1) {
    if (!this->quiet) {
      perror("socket() error");
    }
    return -1;
  }

  // 3. connect
  if (::connect(this->fd, this->servinfo->ai_addr,
                this->servinfo->ai_addrlen) == -1) {
    if (!this->quiet) {
      perror("connect() errror");
    }
    return -1;
  }

  usable = true;
  return 0;
}

int Socket::timedConnect(const std::string& host, int port, double seconds) {
  if (host.empty()) return -1;
  if (port <= 0) return -1;

  // 0. close
  this->close();

  // 1. Get IP
  int status;
  struct addrinfo hints;
  this->servinfo = NULL;            // will point to the results
  memset(&hints, 0, sizeof hints);  // make sure the struct is empty
  hints.ai_family = AF_UNSPEC;      // don't care IPv4 or IPv6
  hints.ai_socktype = SOCK_STREAM;  // TCP stream sockets
  hints.ai_flags = AI_PASSIVE;      // fill in my IP for me
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
  this->fd = ::socket(this->servinfo->ai_family, this->servinfo->ai_socktype,
                      this->servinfo->ai_protocol);
  if (this->fd == -1) {
    if (!this->quiet) {
      perror("socket() error");
    }
    return -1;
  }

  // 3.1 Set non-blocking
  // int socketFlag =fcntl(this->fd, F_GETFL, NULL);
  // if( socketFlag < 0) {
  //   fprintf(stderr, "fcntl F_GETFL error: %s\n", strerror(errno));
  //   return -1;
  // }
  // socketFlag |= O_NONBLOCK;
  if (fcntl(this->fd, F_SETFL, O_NONBLOCK) < 0) {
    fprintf(stderr, "fcntl F_SETFL error: %s\n", strerror(errno));
    return -1;
  }

  // 3.2 Trying to connect with timeout
  status =
      ::connect(this->fd, this->servinfo->ai_addr, this->servinfo->ai_addrlen);
  if ((status == -1) && (errno != EINPROGRESS)) {
    fprintf(stderr, "connect() error: %s\n", strerror(errno));
    close();
    return -1;
  }
  if (status == 0) {
    // connection has succeeded immediately
    usable = true;
    return 0;
  }

  fd_set rfds;
  struct timeval tv;
  tv.tv_sec = (int)seconds;
  tv.tv_usec = (int)((seconds - (int)seconds) * 1000);
  FD_ZERO(&rfds);
  FD_SET(this->fd, &rfds);
  status = select(this->fd + 1, NULL, &rfds, NULL, &tv);

  if (status == 1) {
    int so_error;
    socklen_t len = sizeof so_error;

    getsockopt(this->fd, SOL_SOCKET, SO_ERROR, &so_error, &len);

    if (so_error == 0) {
      // printf("%s:%d is open\n", host.c_str(), port);
    } else {  // error
      fprintf(stderr, "getsockopt() error: %s\n", strerror(errno));
      return -1;
    }
  } else if (status < 0) {  // error happens
    fprintf(stderr, "connect() error %d - %s\n", errno, strerror(errno));
    return -1;
  } else {  // timeout
    fprintf(stderr, "connect() timeout %s:%d\n", host.c_str(), port);
    return -1;
  }

  usable = true;
  return 0;
}

void Socket::close() {
  if (this->fd && this->fd != -1) {  // -1: error or not yet openned
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
  if (!usable) return -1;
  int len = msg.size();
  int nSent = 0;
  while (true) {
    nSent += ::send(this->fd, msg.c_str() + nSent, msg.size(), 0);
    if (nSent == len) {
      break;
    }
    if (nSent < 0) {
      if (!this->quiet) {
        perror("send() error");
      }
      return -1;
    }
  }
  return nSent;
}

int Socket::recv(void* buf, int len) {
  if (!usable) {
    return -1;
  }
  int nRecv = ::recv(this->fd, buf, len, 0);
  if (!quiet && nRecv < 0) {
    perror("recv() error");
  }
  return (nRecv);
}

int Socket::timedRecv(void* buf, int len, double seconds) {
  if (!usable) {
    return -1;
  }

  fd_set readfds;
  FD_ZERO(&readfds);
  FD_SET(this->fd, &readfds);

  struct timeval tv;
  tv.tv_sec = (int)seconds;
  tv.tv_usec = int((seconds - (int)seconds) * 1000);

  int ret = select(this->fd + 1, &readfds, NULL, NULL, &tv);
  if (ret == -1) {  // error occured
    if (!quiet) {
      perror("recv() error in select()");
    }
    return -1;
  } else if (ret == 0) {  // timeout
    if (!quiet) {
      fprintf(stderr, "recv() timeout\n");
    }
    return 0;
  } else {
    return this->recv(buf, len);
  }
}
