#include "Http.h"
#include "Socket.h"
#include "TypeConversion.h"

Http::Http(const std::string& url) {
  // url format:
  // scheme://[user:password@]domain:port/path?query_string#fragment_id
  // ref: https://en.wikipedia.org/wiki/Uniform_resource_locator
  socket = NULL;

  const std::string scheme = "http://";
  if (url.find(scheme) != 0) {
    return;
  }
  std::string noSchemeUrl = url.substr(scheme.size());
  size_t sep = noSchemeUrl.find("/");
  if (sep == std::string::npos) {
    this->domain = noSchemeUrl;
    this->path = "/";
  } else {
    this->domain = noSchemeUrl.substr(0, sep);
    this->path = noSchemeUrl.substr(sep);
  }

  sep = this->domain.find(":");
  if (sep == std::string::npos || sep == this->domain.size()) {
    this->port = 80;
  } else {
    this->port = atoi(this->domain.substr(sep + 1));
    this->domain = this->domain.substr(0, sep);
  }
  
  // check if proxy is used
  char* pProxy = getenv("http_proxy");
  if (pProxy)  {
    // proxy used
    this->proxy = pProxy;
    if (proxy.find(scheme) == 0) {
      this->proxy = proxy.substr(scheme.size());
    }
    size_t sepColon = this->proxy.find(":");
    this->proxyPort = 80;
    if (sepColon != std::string::npos) {
      if (this->proxy.size() == sepColon) { // colon is the last char
        this->proxyPort = 80;
      } else {
        this->proxyPort = atoi(this->proxy.substr(sepColon+1));
      }
      this->proxy = this->proxy.substr(0, sepColon);
    }
    socket = new Socket(this->proxy, this->proxyPort);

    this->request = "GET ";
    this->request += url;
    this->request += " HTTP/1.1\r\nHost: ";
    this->request += this->domain;
    this->request += "\r\n\r\n";
  } else {
    // no proxy
    socket = new Socket(this->domain, port);

    this->request = "GET ";
    this->request += path;
    this->request += "\r\nHost: ";
    this->request += this->domain;
    this->request += "\r\n\r\n";
  }
}

Http::~Http() {
  if (this->socket) {
    delete socket;
    socket = NULL;
  }
}

int Http::read(std::vector<std::string>* content) {
  if (!socket || !socket->isUsable()) {
    return -1;
  }

  socket->send(request);
  int ret = 0;
  content->clear();
  std::string s;
  while (true) {
    ret = socket->timedRecv(buf, bufSize, 2.0);
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
