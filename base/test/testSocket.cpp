#include "Socket.h"
#include <cassert>

int main() {

  {
    Socket s("www.google.com", 80);
    assert(s.isUsable());

    assert(s.send("GET /index.htm HTTP/1.1\r\nHost: google.com\r\n\r\n") > 0);
    char buff[2048];
    int len = s.recv(buff, 2048);

    assert(len > 0);
  }

  {
    Socket s("zhanxw.com", 12345);
    assert(s.send("AA") != 0);
    
  }

  {
    Socket s("unreachable.xxx", 12345);
    assert(s.send("AA") != 0);
    
  }

  return 0;
}
