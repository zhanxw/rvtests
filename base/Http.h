#ifndef HTTP_H
#define HTTP_H

#include <string>
#include <vector>

class Socket;

class HttpResponse {
 public:
  HttpResponse();
  void addLine(const std::string& line);
  void clear();
  size_t size();
  const std::vector<std::string>& getHeader() const { return this->header_; }
  const std::vector<std::string>& getBody() const { return this->body_; }

 private:
  std::vector<std::string> header_;
  std::vector<std::string> body_;
  int status_;
};

class Http {
 public:
  // Make a connect to @param url
  // @param url must be in the format of:
  // http://domain[:port][/path]
  Http(const std::string& url, double timeoutSeconds = -1);
  virtual ~Http();
  /**
   * Read everything and store on lines of contents to @param content
   * If @param maxBodyLines is positive, this functions returns when @param
   * content reaches @param maxBodyLines lines.
   * @return lines of contents read from HTTP, or -1 if error occurs
   */
  int read(HttpResponse* response, int maxBodyLines = -1);
  void enableQuiet();
  void disableQuiet();

 private:
  void stripHeader(std::vector<std::string>* all) const;
#if 0
  bool hasHeader(const std::vector<std::string>& response) const;
  int getStatusCode(const std::vector<std::string>& response) const;
#endif
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
