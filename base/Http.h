#ifndef HTTP_H
#define HTTP_H

#include <string>
#include <vector>

class Socket;
class Http {
 public:
  // Make a connect to @param url
  // @param url must be in the format of:
  // http://domain[:port][/path]
  Http(const std::string& url, double timeoutSeconds = -1);
  virtual ~Http();
  /**
   * Read everything and store on lines of contents to @param content
   * @return lines of contents read from HTTP, or -1 if error occurs
   */
  int read(std::vector<std::string>* content);
  void enableQuiet();
  void disableQuiet();

 private:
  void stripHeader(std::vector<std::string>* all) const;
  bool hasHeader(const std::vector<std::string>& response) const;
  int getStatusCode(const std::vector<std::string>& response) const;

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
