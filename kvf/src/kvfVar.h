// 
// kvfVar.h
// 
// Var pointer class
// 
// David Burgess
// June 1999, June 2004, september 2006
// 

#ifndef KVF_VAR_H
#define KVF_VAR_H

#include "kvfException.h"

namespace KVF {

  // ============================================================ Var Exceptions

  class KVFVarIsNil:public KVFException {
  public:
    KVFVarIsNil(const string & err_str = "", const string fn_str = "")
  :  KVFException(err_str, fn_str) {
      _type_str += "KVFVarIsNil::";
  }};

  class KVFVarSeqIndexError:public KVFException {
  public:
    KVFVarSeqIndexError(const string & err_str = "", const string fn_str = "")
  :  KVFException(err_str, fn_str) {
      _type_str += "KVFVarSeqIndexError::";
  }};
  // -----------------------------------------------------------------------------
  // Class VarTarget
  // -----------------------------------------------------------------------------
  class VarTarget {
    unsigned int ref_count;
  private:
      VarTarget & operator=(const VarTarget & v) {
      return *this;
  } protected:
      VarTarget():ref_count(0) {;
    }                           // default constructor
  public:
    unsigned int get_ref_count() const {
      return ref_count;
    } unsigned int inc_ref_count() {
      return ++ref_count;
    }
    unsigned int dec_ref_count() {
      return --ref_count;
    }
  };

  template < class T > class Var;

  // -----------------------------------------------------------------------------
  // Class Var<T>
  // -----------------------------------------------------------------------------
  template < class T > class Var {
  protected:
    T * p;
  public:

    // copy constructor needed so that ref_count gets incremented
  Var(const Var < T > &var):p(var.p) {
      if (p != 0)
        p->inc_ref_count();
    }

    // constructor from pointer
  Var(T * p_ = 0):p(p_) {
      if (p != 0)
        p->inc_ref_count();
    }

    ~Var() {
      if (p != 0 && p->dec_ref_count() == 0) {
        delete p;
        p = 0;
      }
    }

    T *operator->() const {
      if (p == 0)
        throw KVFVarIsNil();
        return p;
    } Var < T > &operator =(const Var < T > &x) {
      // cout << "Var<T>& operator = ( const Var<T> &x )" << endl;
      if (x.p != 0)
        x.p->inc_ref_count();   // protect against "var = var"
      if (p != 0 && p->dec_ref_count() == 0)
        delete p;
      p = x.p;
      return *this;
    }

    bool is_nil(void) const {
      return (p == 0) ? true : false;
    } bool operator==(const Var < T > &v) const {
      return (p == v.p);
    } bool operator!=(const Var < T > &v) const {
      return (p != v.p);
    }
    // static because does not need access to members of a particular// instance of this class// if T:T2 i.e. T is a derived class of T2 returns a Var<T>// otherwise nil template < class T2 >
    static Var < T > narrow(Var < T2 > object) {
      return Var < T > (dynamic_cast < T * >(object.ptr_()));
    }

    // return Var<T3> from a Var<T> when T3 is a base class of T 
    template < class T3 > Var < T3 > broaden() {
      return p;                 // p is type T*, so initialises via Var<T3>::Var<T3>(T3*)
    }

    template < class T3 > operator  Var < T3 > () {
      return broaden < T3 > ();
    }

    // dereferencing operator

    T & operator*() {
      if (p == 0)
        throw KVFVarIsNil("_var type is nil", "Var::op*()");
      return *p;
    }

    const T & operator*() const {
      if (p == 0)
        throw KVFVarIsNil("_var type is nil", "Var::op*()const");
        return *p;
    } T *ptr_() {
      return p;
    }

    const T *ptr_() const {
      return p;
  }};

  // template<class T>
  // Vostream& operator<<(Vostream& v_os,const Var<T>& var)
  // { const T* ptr=var.ptr_(); v_os << ptr; return v_os;}
  // template<class T>
  // Vistream& operator>>(Vistream& v_is,Var<T>& var)
  // { T* ptr; v_is >> ptr; var=ptr; return v_is;}

  // Var class for container classes with an element index operator
  // Ts is the container class
  // T is the type contained
  // EG: Ts could be Sequence<QString> and T would be QString
  // NOTE: in this implementation T is found from
  // typedef typename Ts::value_type T;
  // which assumes that Ts has vector<T> as a base class, or
  // makes the appropriate definition of value_type


template < class Ts > class SeqVar:public Var < Ts > {

    typedef typename Ts::value_type T;

  public:
    // (default) constructor from pointer
  SeqVar(Ts * p_ = 0):Var < Ts > (p_) {;
    }

    // copy constructor
    SeqVar(const SeqVar < Ts > &svar):Var < Ts > (svar) {;
    }

    // constructor from Var
    // the body of the constructor should really ensure that type Ts
    // inherits from vector<> (or something like it )
    // 
    // hummm not right? SeqVar(const Var<Ts>& var) : Var<Ts>(var) { ; }

    // narrow is static because does not need access to members of a particular
    // instance of this class
    // if Ts:T2 i.e. Ts is a derived class of T2 returns a SeqVar<Ts>
    // otherwise nil

    template < class T2 > static SeqVar < Ts > narrow(Var < T2 > object) {
      return SeqVar < Ts > (dynamic_cast < Ts * >(object.ptr_()));
    }


    T & operator[](int i) {
      if (Var < Ts >::is_nil())
        throw KVFVarIsNil();
      if (i < 0 || i >= Var < Ts >::p->size())
        throw KVFVarSeqIndexError();
      return (*Var < Ts >::p)[i];
    }

    const T & operator[] (int i) const {
      if (Var < Ts >::is_nil())
        throw KVFVarIsNil();
      if (i < 0 || i >= Var < Ts >::p->size())
        throw KVFVarSeqIndexError();
        return (*Var < Ts >::p)[i];
  }};


}                               // end namespace KVF

#endif // #ifndef KVF_VAR_H
