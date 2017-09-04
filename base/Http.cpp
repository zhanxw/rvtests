#include "Http.h"
#include "Socket.h"
#include "TypeConversion.h"
#include "Utils.h"

Http::Http(const std::string& url, double timeoutSeconds) {
  // url format:
  // scheme://[user:password@]domain:port/path?query_string#fragment_id
  // ref: https://en.wikipedia.org/wiki/Uniform_resource_locator
  socket = NULL;
  quiet = false;

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
    // proxy used
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
      // socket shutdown or timeout
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

  // check header
  // NOTE: HTTP header is optional
  if (hasHeader(*content)) {
    // get response code
    if (getStatusCode(*content) != 200) {
      return -1;
    }

    // strip out head lines if any
    stripHeader(content);
  }

  return ((int)content->size());
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
