#ifndef INCLUDE_UTIL_EXCEPTION
#define INCLUDE_UTIL_EXCEPTION

#include <string>

namespace util {

  class Exception {
  public:
    Exception(const std::string& msg = "e") { _msg = msg; }
    virtual ~Exception() {}
    const std::string& GetMsg() const { return _msg; }
  private:
    std::string _msg;
  };

  class FileException : public Exception {
  public:
    FileException(const std::string& msg = "fe") : Exception(msg) {}
  };

  class OutOfMemoryException : public Exception {
  public:
    OutOfMemoryException(const std::string& msg = "oom") : Exception(msg) {}
  };

  class UserReqException : public Exception {
  public:
    UserReqException(const std::string& msg = "ur") : Exception(msg) {}
  };

}

#endif
