#include "Http.h"
#include "Socket.h"
Http::Http(const std::string& url, int port) {
  // url format:
  // scheme://[user:password@]domain:port/path?query_string#fragment_id
  // ref: https://en.wikipedia.org/wiki/Uniform_resource_locator
  socket = NULL;

  const std::string scheme = "http://";
  if (url.find(scheme) == 0) {
    this->domain = url.substr(scheme.size());
  }

  size_t sep = this->domain.find("/");
  if (sep == std::string::npos) {
    this->path = "/";
    this->domain = url;
  } else {
    this->path = this->domain.substr(sep);
    this->domain = this->domain.substr(0, sep);
  }

  socket = new Socket(this->domain, port);
}

Http::~Http() {
  if (this->socket) {
    delete socket;
    socket = NULL;
  }
}

int Http::read(std::vector<std::string>* content) {
  if (!socket->isUsable()) {
    return -1;
  }

  std::string req = "GET ";
  req += path;
  req += "\r\nHost: ";
  req += this->domain;
  req += "\r\n\r\n";

  socket->send(req);
  int ret = 0;
  content->clear();
  std::string s;
  while (true) {
    ret = socket->recv(buf, bufSize);
    if (ret < 0) {
      // error happens
      break;
    }
    if (ret == 0) {
      // socket shutdown
      break;
    }
    for (int i = 0; i < ret; ++i) {
      if (buf[i] == '\r') continue;
      if (buf[i] != '\n') {
        s.push_back(buf[i]);
      } else {
        content->push_back(s);
        s.clear();
      }
    }
  }
  if (!s.empty()) {
    content->push_back(s);
  }
  return ((int)content->size());
}
