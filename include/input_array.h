// input_array.h
// A sample user-defined data type for illustrating ConfigFile
// Operators << and >> are defined to allow writing to and reading from files
// Richard J. Wagner 24 May 2004
// Modified P. Henri 8 June 2011
// corrected by Markidis

#include <iostream>

struct array_int {
  int a, b, c, d, e, f;

    array_int() {
  } array_int(int u1, int u2, int u3, int u4, int u5, int u6):a(u1), b(u2), c(u3), d(u4), e(u5), f(u6) {
  }
  array_int(const array_int & orig):a(orig.a), b(orig.b), c(orig.c), d(orig.d), e(orig.e), f(orig.f) {
  }

  array_int & operator=(const array_int & orig) {
    a = orig.a;
    b = orig.b;
    c = orig.c;
    d = orig.d;
    e = orig.e;
    f = orig.f;
    return *this;
  }
};


inline std::ostream & operator<<(std::ostream & os, const array_int & t) {

  os << t.a << " " << t.b << " " << t.c << " " << t.d << " " << t.e << " " << t.f;
  return os;
}


inline std::istream & operator>>(std::istream & is, array_int & t) {

  is >> t.a >> t.b >> t.c >> t.d >> t.e >> t.f;
  return is;
}


struct array_double {
  double a, b, c, d, e, f;

  array_double() {
  } array_double(int u1, int u2, int u3, int u4, int u5, int u6):a(u1), b(u2), c(u3), d(u4), e(u5), f(u6) {
  }
  array_double(const array_double & orig):a(orig.a), b(orig.b), c(orig.c), d(orig.d), e(orig.e), f(orig.f) {
  }

  array_double & operator=(const array_double & orig) {
    a = orig.a;
    b = orig.b;
    c = orig.c;
    d = orig.d;
    e = orig.e;
    f = orig.f;
    return *this;
  }
};


inline std::ostream & operator<<(std::ostream & os, const array_double & t) {

  os << t.a << " " << t.b << " " << t.c << " " << t.d << " " << t.e << " " << t.f;
  return os;
}


inline std::istream & operator>>(std::istream & is, array_double & t) {
  is >> t.a >> t.b >> t.c >> t.d >> t.e >> t.f;
  return is;
}


struct array_bool {
  bool a, b, c, d, e, f;

  array_bool() {
  } array_bool(int u1, int u2, int u3, int u4, int u5, int u6):a(u1), b(u2), c(u3), d(u4), e(u5), f(u6) {
  }
  array_bool(const array_bool & orig):a(orig.a), b(orig.b), c(orig.c), d(orig.d), e(orig.e), f(orig.f) {
  }

  array_bool & operator=(const array_bool & orig) {
    a = orig.a;
    b = orig.b;
    c = orig.c;
    d = orig.d;
    e = orig.e;
    f = orig.f;
    return *this;
  }
};


inline std::ostream & operator<<(std::ostream & os, const array_bool & t) {
  os << t.a << " " << t.b << " " << t.c << " " << t.d << " " << t.e << " " << t.f;
  return os;
}


inline std::istream & operator>>(std::istream & is, array_bool & t) {
  is >> t.a >> t.b >> t.c >> t.d >> t.e >> t.f;
  return is;
}
