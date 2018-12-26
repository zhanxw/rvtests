#include "Http.h"

#include <stdlib.h>  //getenv
#include "Socket.h"
#include "TypeConversion.h"
#include "Utils.h"

HttpResponse::HttpResponse() {
  // status_ table
  // 0 : uninitialized
  // 1 : header
  // 2 : body
  status_ = 0;
}

void HttpResponse::addLine(const std::string& line) {
  if (status_ == 0) {
    if (line.size()) {
      if (line.substr(0, 5) == "HTTP/") {
        status_ = 1;
        header_.push_back(line);
      } else {
        status_ = 2;
        body_.push_back(line);
      }
    }
  } else if (status_ == 1) {
    if (line.size()) {
      header_.push_back(line);
    } else {
      status_ = 2;
    }
  } else if (status_ == 2) {
    body_.push_back(line);
  }
}

void HttpResponse::clear() {
  header_.clear();
  body_.clear();
  status_ = 0;
}

size_t HttpResponse::size() { return (header_.size() + body_.size()); }

//////////////////////////////////////////////////

Http::Http(const std::string& url, double timeoutSeconds) {
  // url format:
  // scheme://[user:password@]domain:port/path?query_string#fragment_id
  // ref: https://en.wikipedia.org/wiki/Uniform_resource_locator
  socket = NULL;
  quiet = true;
  if (std::getenv("HTTP_DEBUG")) {
    quiet = false;
  }

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
  if (pProxy) {
    if (!this->quiet) {
      fprintf(stderr, "Use proxy [ %s ]\n", pProxy);
    }
    this->proxy = pProxy;
    if (proxy.find(scheme) == 0) {
      this->proxy = proxy.substr(scheme.size());
    }
    size_t sepColon = this->proxy.find(":");
    this->proxyPort = 80;
    if (sepColon != std::string::npos) {
      if (this->proxy.size() == sepColon) {  // colon is the last char
        this->proxyPort = 80;
      } else {
        this->proxyPort = atoi(this->proxy.substr(sepColon + 1));
      }
      this->proxy = this->proxy.substr(0, sepColon);
    }
    if (timeoutSeconds > 0) {
      socket = new Socket(this->proxy, this->proxyPort, timeoutSeconds);
    } else {
      socket = new Socket(this->proxy, this->proxyPort);
    }

    this->request = "GET ";
    this->request += url;
    this->request += " HTTP/1.1\r\nHost: ";
    this->request += this->domain;
    this->request += "\r\n\r\n";
  } else {
    // no proxy
    if (!this->quiet) {
      fprintf(stderr, "No proxy used\n");
    }
    if (timeoutSeconds > 0) {
      socket = new Socket(this->domain, port, timeoutSeconds);
    } else {
      socket = new Socket(this->domain, port);
    }

    this->request = "GET ";
    this->request += path;
    this->request += "\r\nHost: ";
    this->request += this->domain;
    this->request += "\r\n\r\n";
  }
  if (this->quiet) {
    this->enableQuiet();
  } else {
    this->disableQuiet();
  }
}

Http::~Http() {
  if (this->socket) {
    delete socket;
    socket = NULL;
  }
}

int Http::read(HttpResponse* response, int maxBodyLines) {
  if (!socket || !socket->isUsable()) {
    return -1;
  }

  socket->send(request);
  int ret = 0;
  std::string s;
  response->clear();
  while (true) {
    ret = socket->timedRecv(buf, bufSize, 2.0);
    if (ret < 0) {
      // error happens
      if (!this->quiet) {
        fprintf(stderr, "timedRecv failed with ret = [ %d ]\n", ret);
      }
      break;
    }
    if (ret == 0) {
      // socket shutdown or timeout
      if (!this->quiet) {
        fprintf(
            stderr,
            "timedRecv receives nothing as socket shuts down or is timedout\n");
      }
      break;
    }
    for (int i = 0; i < ret; ++i) {
      if (buf[i] == '\r') continue;
      if (buf[i] != '\n') {
        s.push_back(buf[i]);
      } else {
        response->addLine(s);
        s.clear();
      }
    }
    if (!this->quiet) {
      fprintf(stderr, "timedRecv receives [ %d ] bytes: [ ", ret);
      for (int i = 0; i < ret; ++i) {
        if (buf[i] == '\r') {
          fprintf(stderr, "\\r");
        } else if (buf[i] == '\n') {
          fprintf(stderr, "\\n");
        } else {
          fprintf(stderr, "%c", buf[i]);
        }
      }
      fprintf(stderr, " ]\n");
    }

    if (maxBodyLines > 0 &&
        response->getBody().size() >= (size_t)maxBodyLines) {
      break;
    }
  }
  if (!s.empty()) {  // remaining bytes read, so need to place it.
    response->addLine(s);
  }

  return ((int)response->size());
}

void Http::enableQuiet() {
  if (this->socket) {
    this->socket->enableQuiet();
  }
}

void Http::disableQuiet() {
  if (this->socket) {
    this->socket->disableQuiet();
  }
}

void Http::stripHeader(std::vector<std::string>* all) const {
  std::vector<std::string>& content = *all;
  size_t sep = 0;
  for (size_t i = 1; i != content.size(); ++i) {
    if (!content[i].size()) {
      sep = i;
      break;
    }
  }
  if (sep > 0) {
    size_t n = content.size() - sep - 1;
    for (size_t i = 0; i != n; ++i) {
      (content)[i] = content[i + sep + 1];
    }
    content.resize(n);
  }
}

#if 0
bool Http::hasHeader(const std::vector<std::string>& response) const {
  if (response.size() && response[0].substr(0, 5) == "HTTP/") {
    return true;
  }
  return false;
}

int Http::getStatusCode(const std::vector<std::string>& response) const {
  const std::string& line = response[0];
  std::vector<std::string> res;
  if (stringNaturalTokenize(line, " \t", &res) < 2) {
    return -1;
  }
  int code = -1;
  if (str2int(res[1], &code)) return code;
  return -1;
}
#endif
