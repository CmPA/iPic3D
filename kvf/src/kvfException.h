// 
// KVFException.h
// 
// Exception class
// 
// David Burgess
// September 2006
// 

/* ! \file Header file for KVFException class.
 * 
 */

#ifndef KVF_EXCEPTION_H
#define KVF_EXCEPTION_H

#include <sstream>
#include <iostream>
#include <string>
#include <vector>

namespace KVF {

  using namespace std;

  class KVFException {
  protected:
    bool _can_recover;          // /< Flag indicating recovery possible
    string _err_str;            // /< Exception error message string
    string _fn_str;             // /< Function name throwing exception
    string _type_str;           // /< Exception type string, including inheritance
      vector < string > _err_msgs;  // /< list of associated messages
    int _sys_errno;             // /< system error number from bad system call
    string _sys_err_str;        // /< system error message from bad system call
  public:
      KVFException(void):_type_str("KVFException") {;
    } KVFException(const string & err_str, const string fn_str = "", int sys_errno = 0)
  :  _can_recover(true), _err_str(err_str), _fn_str(fn_str), _type_str("MdException::"), _sys_errno(sys_errno) {
      if (sys_errno != 0)
        _sys_err_str = strerror(sys_errno);
    } KVFException(const KVFException & e):_can_recover(e._can_recover), _err_str(e._err_str), _fn_str(e._fn_str), _type_str(e._type_str), _err_msgs(e._err_msgs), _sys_errno(e._sys_errno), _sys_err_str(e._sys_err_str) {;
    }

    void append_err_str(const string & str) {
      _err_str += str;
    }
    void prepend_err_str(const string & str) {
      _err_str = str + _err_str;
    }
    void push_err_msg(const string & err_msg) {
      _err_msgs.push_back(err_msg);
    }
    void diag_cout() {
      cout << _type_str << " " << _err_str << " in function " << _fn_str << "\n";
      if (_sys_errno != 0)
        cout << "sys errno: " << _sys_errno << " [" << _sys_err_str << "]" << "\n";
      int n_msgs = _err_msgs.size();
      for (int i = 0; i < n_msgs; ++i)
        cout << ".. " << _err_msgs[i] << "\n";
    }

  };

}                               // end namespace KVF

#endif // #ifndef KVF_EXCEPTION_H
