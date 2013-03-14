// 
// Naming.h
// 
// KVF - Name Context
// 
// David Burgess
// June 1999, June 2004, September 2006
// 

#ifndef KVF_NAMING_H_
#define KVF_NAMING_H_

#include <regex.h>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "kvfException.h"
#include "kvfVar.h"

namespace KVF {

  using namespace std;

  // -------------------------------------------------------------------------

  class NameComponent {
    string identifier;
    string kind_attrib;

    friend class Name;

  public:
      NameComponent(string id, string kind = "")
  :  identifier(id), kind_attrib(kind) {;
    } NameComponent(const NameComponent & a)
    : identifier(a.identifier), kind_attrib(a.kind_attrib) {;
    } string get_id(void) const {
      return identifier;
    } string set_id(const string & id) {
      identifier = id;
      return identifier;
    }
    bool is_valid_id(void) const {
      return identifier != "";
    } bool operator==(const string & a) const {
      return a == identifier;
  }};

  class Name {
    vector < NameComponent > components;
  public:
    Name(const string & name, string separators = "/");
    Name(const Name & name);
      Name(const NameComponent & name_component);

    int num_components(void) const {
      return components.size();
    } NameComponent & get_component(int i);
    const NameComponent & get_component(int i) const;
    void delete_component(int i);

    const NameComponent & operator[] (int i) const;

      NameComponent & operator[] (int i);

    string srep(const string & separator = "/") const;
    string strrep(const string & separator = "/") const {
      return srep();
    } string s_rep(const string & separator = "/") const {
      return srep();
    } void diag_print() const {
      cout << "KVF::Name diag_print\n";
      for (int i = 0; i < components.size(); ++i) {
        cout << "[" << i << "] <" << components[i].identifier << ">\n";
  }}};

  // ===================================================== EXCEPTION CLASSES

  class KVFNamingException:public KVFException {
  public:
    KVFNamingException(const string & err_str = "", const string fn_str = "")
  :  KVFException(err_str, fn_str) {
      _type_str += "KVFNamingException::";
  }};

  class KVFNameInvalid:public KVFNamingException {
  public:
    KVFNameInvalid(const string & err_str = "", const string fn_str = "")
  :  KVFNamingException(err_str, fn_str) {
      _type_str += "KVFNameInvalid:";
  }};

  class KVFNamingError:public KVFNamingException {
  public:
    KVFNamingError(const string & err_str = "", const string fn_str = "")
  :  KVFNamingException(err_str, fn_str) {
      _type_str += "KVFNamingError:";
  }};

  class KVFNameAlreadyBound:public KVFNamingException {
  public:
    KVFNameAlreadyBound(const string & err_str = "", const string fn_str = "")
  :  KVFNamingException(err_str, fn_str) {
      _type_str += "KVFNameAlreadyBound:";
    };
  };

  class KVFNameNotFound:public KVFNamingException {
    Name _not_found_name;
  public:
      KVFNameNotFound(const Name & n)
    : KVFNamingException("NameNotFound:"), _not_found_name(n) {
      _type_str += "KVFNameNotFound:";
    } Name name() {
      return _not_found_name;
    }
  };


  // ========================================================================



  typedef enum BindingTypeE { N_OBJECT, N_CONTEXT } BindingType;

  class Binding {
  public:
    Name name;
    BindingType type;
    string desc;
      Binding(void):name(""), type(N_OBJECT), desc("") {;
    } Binding(const Binding & a):name(a.name), type(a.type), desc(a.desc) {;
    }
    Binding(Name n, BindingType t, string d)
  :  name(n), type(t), desc(d) {;
    }
  };

  typedef vector < Binding > BindingList;

  template < class T > class NamingContext;

template < class T > class NamingContext:public VarTarget {
  public:
    typedef Var < NamingContext > _Var;
    typedef NamingContext *_Ptr;
    typedef const NamingContext *_CPtr;

  private:

    // exception which uses _Var, so has to be in template class
    // .. make it public?
    // NB: should use _var instead of ptr, but ownership problems,
    // since initialized from a "this" -- tobesolved?

  class NameNotFoundT:public KVFNameNotFound {
      NamingContext *_naming_context_ptr;
    public:
      NameNotFoundT(NamingContext * nc, Name n)
    :  _naming_context_ptr(nc), KVFNameNotFound(n) {;
      }

      NamingContext *naming_context(void) {
        return _naming_context_ptr;
      }
    };

    class BindingEntry {
    public:
      T _obj;
      _Var _nc;
      BindingType _type;
      string _desc;

        BindingEntry(void):_obj(), _nc(0), _type(N_OBJECT), _desc("") {;
      } BindingEntry(const BindingEntry & a):_type(a._type), _desc(a._desc) {
        if (_type == N_CONTEXT)
          _nc = a._nc;
        else
          _obj = a._obj;
      }

      BindingEntry(const T & obj, const string & desc)
      : _obj(obj), _nc(0), _type(N_OBJECT), _desc(desc) {;
      }

      BindingEntry(_Var nc, const string & desc)
      : _nc(nc), _type(N_CONTEXT), _desc(desc) {;
      }

      const BindingEntry & operator=(const BindingEntry & a) {
        _type = a._type;
        _desc = a._desc;
        if (_type == N_CONTEXT)
          _nc = a._nc;
        else
          _obj = a._obj;
        return *this;
      }

      ~BindingEntry(void) {;
      }

      T & n_object() {
        return _obj;
      }
    };

    typedef map < string, BindingEntry > BindingMap;

    BindingMap bindings;

  public:

    NamingContext(void) {;
    }
    // NamingContext( const NameContext& a );

    ~NamingContext(void) {;
    }
    // { cout << "~NamingContext"<<endl; }

    class recursive_iterator:public BindingMap::iterator, public VarTarget {
    protected:
      typedef Var < recursive_iterator > iterator_var;
      iterator_var _parent;
      _Ptr _nc;
    public:
        recursive_iterator():BindingMap::iterator(), _nc(0), _parent(0) {;
      } recursive_iterator(_Ptr nc, const typename BindingMap::iterator & itr):_nc(nc), BindingMap::iterator(itr), _parent(0) {;
      }
      recursive_iterator(const recursive_iterator & itr)
      : _nc(itr._nc), _parent(itr._parent), BindingMap::iterator(itr) {;
      }
      ~recursive_iterator() {;
      }

      recursive_iterator & operator=(const recursive_iterator & itr) {
        _nc = itr._nc;
        _parent = itr._parent;
        BindingMap::iterator::operator=(itr);
        return *this;
      }
      recursive_iterator & operator++() { // if not past-the-end, get entry type
        if ((*this) != _nc->recursive_end()) {
          iterator itr(_nc, *this);
          if (itr.type() == N_CONTEXT) {
            // if N_CONTEXT get beginnining of next naming context
            recursive_iterator tmp = *this;
            (*this) = (itr.naming_context())->recursive_begin();
            // store details of originating point
            _parent = new recursive_iterator(tmp);
          }
          else
            BindingMap::iterator::operator++(); // simple increment if N_OBJECT
        }
        // if reach end of Naming Context go up a level if possible
        while ((*this) == _nc->recursive_end() && !_parent.is_nil()) {
          recursive_iterator tmp = *_parent;
          (*this) = tmp;
          BindingMap::iterator::operator++();
        }
        return *this;
      }
      recursive_iterator operator++(int) {
        recursive_iterator tmp = *this;
        ++*this;
        return tmp;
      }
      recursive_iterator & operator--() { // if at beginning of a BindingMap, get previous entry off
        if ((*this) == _nc->recursive_begin()) {
          if (!_parent.is_nil()) {
            recursive_iterator tmp = *_parent;
            (*this) = tmp;
          }

        }
        else {                  // get previous entry on bindingmap
          BindingMap::iterator::operator--();
          // get end of bindingmap pointed to by Naming Context
          iterator itr(_nc, *this);
          if (itr.type() == N_CONTEXT) {
            recursive_iterator tmp = *this;
            (*this) = (itr.naming_context())->recursive_end();
            _parent = new recursive_iterator(tmp);
            --*this;
          }
        }
        return *this;
      }

      recursive_iterator operator--(int) {
        recursive_iterator tmp = *this;
        --*this;
        return tmp;
      }

      BindingType type() const {
        // Return N_CONTEXT or N_OBJECT if not past-the-end 
        if ((*this) == _nc->recursive_end())
          throw KVFNamingError();
          return (*this)->second._type;
      } string name(string separator = "/") const {
        // get full Naming Context name
        string binding_name;
        if (!_parent.is_nil())
            binding_name = _parent->name(separator) + separator;
        if ((*this) != _nc->recursive_end())
            binding_name += (*this)->first;
          return binding_name;
      } _Ptr nc() const {
        return _nc;
    } private:
      void rdepth(int &d) const {
        if (!_parent.is_nil()) {
          d++;
          _parent->rdepth(d);
    }} public:
      int depth() const {
        int d = 0;
        if (!_parent.is_nil())
            rdepth(d);
          return d;
      } bool match(const string & pattern) const {
        // check whether name matches pattern
        if ((*this) == _nc->recursive_end())
          return true;
        if (!pattern.size())
          return true;
        regex_t re;
        if (regcomp(&re, pattern.c_str(), REG_EXTENDED | REG_NOSUB) != 0)
            return false;
        string str = name();
        int status = regexec(&re, str.c_str(), (size_t) 0, NULL, 0);
          regfree(&re);
        if (status != 0)
            return false;
          return true;
    }};

    class const_recursive_iterator:public BindingMap::const_iterator, public VarTarget {
    protected:
      typedef Var < const_recursive_iterator > const_iterator_var;
      const_iterator_var _parent;
      _CPtr _nc;
    public:
        const_recursive_iterator()
      : BindingMap::const_iterator(), _nc(0), _parent(0) {;
      } const_recursive_iterator(_CPtr nc, const typename BindingMap::const_iterator & citr)
      : _nc(nc), BindingMap::const_iterator(citr), _parent(0) {;
      }
      const_recursive_iterator(const const_recursive_iterator & citr)
      : _nc(citr._nc), _parent(citr._parent), BindingMap::const_iterator(citr) {;
      }
      ~const_recursive_iterator() {;
      }

      const_recursive_iterator & operator=(const const_recursive_iterator & citr) {
        _parent = citr._parent;
        _nc = citr._nc;
        BindingMap::const_iterator::operator=(citr);
        return *this;
      }
      const_recursive_iterator & operator++() { // if not past-the-end, get entry type
        if ((*this) != _nc->recursive_end()) {
          const_iterator citr(_nc, *this);
          if (citr.type() == N_CONTEXT) {
            // if N_CONTEXT get beginnining of next naming context
            const_recursive_iterator tmp = *this;
            const NamingContext & nc = *citr.naming_context();
            (*this) = nc.recursive_begin();
            // store details of originating point
            _parent = new const_recursive_iterator(tmp);
          }
          else
            BindingMap::const_iterator::operator++();
        }
        // if reach end of Naming Context go up a level if possible
        while ((*this) == _nc->recursive_end() && !_parent.is_nil()) {
          (*this) = *_parent;
          BindingMap::const_iterator::operator++();
        }

        return *this;
      }

      const_recursive_iterator operator++(int) {
        const_recursive_iterator tmp = *this;
        ++*this;
        return tmp;
      }

      const_recursive_iterator & operator--() {
        // if at beginning of a BindingMap, get previous entry off stack
        if ((*this) == _nc->recursive_begin()) {
          if (!_parent.is_nil())
            (*this) = *_parent;
        }
        else {
          // get previous entry on bindingmap
          BindingMap::const_iterator::operator--();
          // get end of bindingmap pointed to if by Naming Context
          const_iterator citr(_nc, *this);
          if (citr.type() == N_CONTEXT) {
            const_recursive_iterator tmp = *this;
            const NamingContext & nc = *citr.naming_context();
            (*this) = nc.recursive_end();
            _parent = new const_recursive_iterator(tmp);
            --*this;
          }
        }

        return *this;
      }

      const_recursive_iterator operator--(int) {
        const_recursive_iterator tmp = *this;
        --*this;
        return tmp;
      }

      BindingType type() const {
        // Return N_CONTEXT or N_OBJECT if not past-the-end 
        if ((*this) == _nc->recursive_end())
          throw KVFNamingError();
          return (*this)->second._type;
      } string name(string delim = "/") const {
        // get full Naming Context name
        string binding_name;
        if (!_parent.is_nil())
            binding_name = _parent->name() + delim;
        if ((*this) != _nc->recursive_end())
            binding_name += (*this)->first;
          return binding_name;
      } _CPtr nc() const {
        return _nc;
    } private:
      void rdepth(int &d) const {
        if (!_parent.is_nil()) {
          d++;
          _parent->rdepth(d);
    }} public:
      int depth() const {
        int d = 0;
        if (!_parent.is_nil())
            rdepth(d);
          return d;
      } bool match(const string & pattern) const {
        // check whether name matches pattern
        if ((*this) == _nc->recursive_end())
          return true;
        if (!pattern.size())
          return true;
        regex_t re;
        if (regcomp(&re, pattern.c_str(), REG_EXTENDED | REG_NOSUB) != 0)
            return false;
        string str = name();
        int status = regexec(&re, str.c_str(), (size_t) 0, NULL, 0);
          regfree(&re);
        if (status != 0)
            return false;
          return true;
    }};                         // end class const_recursive_iterator

    class const_iterator;

    class iterator:public BindingMap::iterator {
      _Ptr _nc;

      friend class NamingContext < T >::const_iterator;

    public:
        iterator():_nc(0), BindingMap::iterator() {;
      } iterator(_Ptr nc, const typename BindingMap::iterator & itr)
      : _nc(nc), BindingMap::iterator(itr) {;
      }
      ~iterator() {;
      }

      string name() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("iterator::name() off end");
          return (*this)->first;
      } BindingType type() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("iterator::type() off end");
          return (*this)->second._type;
      } T object() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("iterator::object() off end");
        if ((*this)->second._type != N_OBJECT)
          throw KVFNamingError("iterator::object() Not an object");
          return (*this)->second._obj;
      } _Var naming_context() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("iterator::object() off end");
        if ((*this)->second._type != N_CONTEXT)
          throw KVFNamingError("iterator::object() Not an object");
          return (*this)->second._nc;
    }};                         // end class iterator

    class const_iterator:public BindingMap::const_iterator {
      _CPtr _nc;
    public:
        const_iterator():_nc(0), BindingMap::const_iterator() {;
      } const_iterator(_CPtr nc, const typename BindingMap::const_iterator & itr)
      : _nc(nc), BindingMap::const_iterator(itr) {;
      }
      const_iterator(const iterator & itr):_nc(itr._nc), BindingMap::const_iterator(itr) {;
      }
      ~const_iterator() {;
      }

      const_iterator & operator=(const const_iterator & citr) {
        _nc = citr._nc;
        BindingMap::const_iterator::operator=(citr);
        return *this;
      }

      string name() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("recursion off end", "NamingContext::const_iterator::name()");
          return (*this)->first;
      } BindingType type() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("recursion off end", "NamingContext::const_iterator::type()");
          return (*this)->second._type;
      } T object() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("recursion off end", "NamingContext::const_iterator::object()");
        if ((*this)->second._type != N_OBJECT)
          throw KVFNamingError("recursion off end", "NamingContext::const_iterator::object()");
          return (*this)->second._obj;
      } _Var naming_context() const {
        if ((*this) == _nc->end())
          throw KVFNamingError("recursion off end", "NamingContext::const_iterator::naming_context()");
        if ((*this)->second._type != N_CONTEXT)
          throw KVFNamingError("recursion off end", "NamingContext::const_iterator::naming_context()");
          return (*this)->second._nc;
    }};                         // end class const_iterator

    class matched_iterator:public recursive_iterator {
      string _pattern;

    public:
        matched_iterator(const string & pattern = "") {;
    } matched_iterator(_Ptr nc, const typename BindingMap::iterator & itr, const string & pattern = ""):recursive_iterator(nc, itr), _pattern(pattern) {
        while (!recursive_iterator::match(pattern))
          recursive_iterator::operator++();
      }
    matched_iterator(const recursive_iterator & itr, const string & pattern = ""):
      recursive_iterator(itr), _pattern(pattern) {
        while (!recursive_iterator::match(pattern))
          recursive_iterator::operator++();
      }
      matched_iterator(const matched_iterator & mitr):recursive_iterator(mitr), _pattern(mitr._pattern) {;
      }

      matched_iterator & operator=(const recursive_iterator & itr) {
        recursive_iterator::operator=(itr);
        while (!recursive_iterator::match(_pattern))
          recursive_iterator::operator++();
        return *this;
      }

      matched_iterator & operator++() {
        recursive_iterator::operator++();
        while (!recursive_iterator::match(_pattern))
          recursive_iterator::operator++();
        return *this;
      }

      matched_iterator operator++(int) {
        matched_iterator tmp = *this;
        ++*this;
        return tmp;
      }

      /* Not right? matched_iterator& operator--() { if((*this)!=_nc->recursive_begin(_pattern)) recursive_iterator::operator--(); while(!recursive_iterator::match(_pattern)) recursive_iterator::operator--(); return *this; } */

      matched_iterator operator--(int) {
        matched_iterator tmp = *this;
        --*this;
        return tmp;
      }
    };                          // end class matched_iterator

    class const_matched_iterator:public const_recursive_iterator {
      string _pattern;

    public:
        const_matched_iterator(const string & pattern = "") {;
    } const_matched_iterator(_CPtr nc, const typename BindingMap::const_iterator & itr, const string & pattern = ""):const_recursive_iterator(nc, itr), _pattern(pattern) {
        while (!const_recursive_iterator::match(pattern))
          const_recursive_iterator::operator++();
      }
    const_matched_iterator(const const_recursive_iterator & itr, const string & pattern = ""):const_recursive_iterator(itr), _pattern(pattern) {
        while (!const_recursive_iterator::match(pattern))
          const_recursive_iterator::operator++();
      }
      const_matched_iterator(const const_matched_iterator & mitr):const_recursive_iterator(mitr), _pattern(mitr._pattern) {;
      }

      const_matched_iterator & operator=(const const_recursive_iterator & itr) {
        const_recursive_iterator::operator=(itr);
        while (!const_recursive_iterator::match(_pattern))
          const_recursive_iterator::operator++();
        return *this;
      }

      const_matched_iterator & operator++() {
        const_recursive_iterator::operator++();
        while (!const_recursive_iterator::match(_pattern))
          const_recursive_iterator::operator++();
        return *this;
      }

      const_matched_iterator operator++(int) {
        const_matched_iterator tmp = *this;
        ++*this;
        return tmp;
      }

      /* Not right? const_matched_iterator& operator--() { if((*this)!=_nc->recursive_begin(_pattern)) const_recursive_iterator::operator--(); while(!const_recursive_iterator::match(_pattern)) const_recursive_iterator::operator--(); return *this; }
       * 
       */
      const_matched_iterator operator--(int) {
        const_matched_iterator tmp = *this;
        --*this;
        return tmp;
      }
    };                          // class const_matched_iterator

    iterator begin() {
      typename BindingMap::iterator itr = bindings.begin();
      return iterator(this, itr);
    }

    iterator end() {
      typename BindingMap::iterator itr = bindings.end();
      return iterator(this, itr);
    }

    const_iterator begin() const {
      typename BindingMap::const_iterator citr = bindings.begin();
        return const_iterator(this, citr);
    } const_iterator end() const {
      typename BindingMap::const_iterator citr = bindings.end();
        return const_iterator(this, citr);
    } matched_iterator recursive_begin(const string & pattern = "") {
      typename BindingMap::iterator itr = bindings.begin();
      return matched_iterator(this, itr, pattern);
    }

    recursive_iterator recursive_end() {
      typename BindingMap::iterator itr = bindings.end();
      return recursive_iterator(this, itr);
    }

    const_matched_iterator recursive_begin(const string & pattern = "") const {
      typename BindingMap::const_iterator citr = bindings.begin();
        return const_matched_iterator(this, citr, pattern);
    } const_recursive_iterator recursive_end() const {
      typename BindingMap::const_iterator citr = bindings.end();
        return const_recursive_iterator(this, citr);
    } T resolve(const string & n) {
      return resolve(Name(n));
    }

    void bind(const string & n, const T & obj, string desc = "") {
      bind(Name(n), obj, desc);
    }

    void bindx(const string & n, const T & obj, string desc = "") {
      bindx(Name(n), obj, desc);
    }

    void bind_context(const string & n, _Var nc, string desc = "") {
      bind_context(Name(n), nc, desc);
    }

    void unbind(const string & n) {
      unbind(Name(n));
    }

    bool context_is_empty(const string & n) {
      return context_is_empty(Name(n));
    }

    _Var resolve_context(const string & n) {
      return resolve_context(Name(n));
    }

    T resolve(Name n) {
      int n_components = n.num_components();
      if (n_components == 0)
        throw KVFNameInvalid();
      if (n_components == 1) {
        // simple name
        if (!n[0].is_valid_id())
          throw KVFNameInvalid();
        typename BindingMap::const_iterator binding = bindings.find(n[0].get_id());
        if (binding != bindings.end())
          if (binding->second._type == N_OBJECT)
            return binding->second._obj;
          else
            throw KVFNameInvalid();
        else
          throw NameNotFoundT(this, n);
      }
      else {
        // complex name
        Name simple_name = Name(n[n_components - 1]);
        n.delete_component(n_components - 1);
        return resolve_context(n)->resolve(simple_name);
      }
    }                           // end resolve

    _Var resolve_context(Name n) {
      int n_components = n.num_components();
      if (n_components == 0)
        throw KVFNameInvalid();
      if (n_components == 1) {
        // simple name
        if (!n[0].is_valid_id())
          throw KVFNameInvalid();
        typename BindingMap::const_iterator binding = bindings.find(n[0].get_id());
        if (binding != bindings.end()) {
          if (binding->second._type == N_CONTEXT)
            return binding->second._nc;
          else
            throw KVFNameInvalid();
        }
        else
          throw NameNotFoundT(this, n);
      }
      else {
        // complex name
        Name simple_name = Name(n[n_components - 1]);
        n.delete_component(n_components - 1);
        _Var nc = resolve_context(n);
        return nc->resolve_context(simple_name);
      }
    }                           // end resolve_context

    // bindx will do bind_new_context (repeatedly) if not all components
    // of the name have a corresponding context
    void bindx(Name n, const T & obj, string desc = "") {
      bool bind_done = false;
      while (!bind_done) {
        try {
          bind(n, obj, desc);
          bind_done = true;
        }
        catch(NameNotFoundT & not_found) {
          not_found.naming_context()->bind_new_context(not_found.name());
        }
      }
    }

    void bind(Name n, const T & obj, string desc = "") {
      int n_components = n.num_components();
      if (n_components == 0)
        throw KVFNameInvalid();
      if (n_components == 1) {
        // simple name
        typename BindingMap::const_iterator binding = bindings.find(n[0].get_id());
        if (binding != bindings.end())
          throw KVFNameAlreadyBound();
        else {
          BindingEntry binding_entry(obj, desc);
          bindings[n[0].get_id()] = binding_entry;
        }
      }
      else {
        // complex name
        Name simple_name = Name(n[n_components - 1]);
        n.delete_component(n_components - 1);
        resolve_context(n)->bind(simple_name, obj, desc);
      }
    }

    void bind_context(Name n, _Var nc, string desc = "") {
      int n_components = n.num_components();
      if (n_components == 0)
        throw KVFNameInvalid();
      if (n_components == 1) {
        // simple name
        if (!n[0].is_valid_id())
          throw KVFNameInvalid();
        typename BindingMap::const_iterator binding = bindings.find(n[0].get_id());
        if (binding != bindings.end())
          throw KVFNameAlreadyBound();
        else
        bindings[n[0].get_id()] = BindingEntry(nc, desc);
      }
      else {
        // complex name
        Name simple_name = Name(n[n_components - 1]);
        n.delete_component(n_components - 1);
        resolve_context(n)->bind_context(simple_name, nc, desc);
      }
    }

    void rebind(Name n, T obj, string desc = "");
    void rebind_context(Name n, _Var nc, string desc = "");

    void bind_new_context(Name n, string desc = "") {
      _Var new_nc = new NamingContext();
      bind_context(n, new_nc, desc);
    }

    bool context_is_empty(Name n) {
      _Var nc = resolve_context(n);
      return nc->bindings.size() == 0;
    }

    bool is_empty() {
      return bindings.size() == 0;
    }

    void unbind(Name n) {
      int n_components = n.num_components();
      if (n_components == 0)
        throw KVFNameInvalid("Trying to unbind for zero n_components", "NamingContext::unbind");
      if (n_components == 1) {
        // simple name
        typename BindingMap::iterator binding = bindings.find(n[0].get_id());
        if (binding != bindings.end()) {

          bindings.erase(binding);

        }
        else
          throw NameNotFoundT(this, n);
      }
      else {
        // complex name
        Name simple_name = Name(n[n_components - 1]);
        n.delete_component(n_components - 1);
        _Var nc = resolve_context(n);
        nc->unbind(simple_name);
      }
    }                           // end unbind

    void list(BindingList & bl) const {
      bl.clear();
      for (typename BindingMap::const_iterator binding = bindings.begin(); binding != bindings.end(); ++binding)
        bl.push_back(Binding(Name(binding->first), binding->second._type, binding->second._desc));
    } // end list void list_names_r(
                                     vector < string > &sl, string name_root = "", string separator = "/") {
      for (typename BindingMap::const_iterator binding = bindings.begin(); binding != bindings.end(); ++binding) {

        string binding_name = name_root + binding->first;
        if (binding->second._type == N_OBJECT)
          sl.push_back(binding_name);
        else if (binding->second._type == N_CONTEXT)
          binding->second._nc->list_names_r(sl, binding_name + separator, separator);
      }
    }
    void diag_print() {
      recursive_iterator itr;
      for (itr = recursive_begin(); itr != recursive_end(); itr++) {
        cout << itr.name();
        if (itr.type() == N_CONTEXT)
          cout << "/";
        cout << "  [depth: " << itr.depth() << "]" << endl;
      }
    }

  };

}                               // end namespace KVF

# endif // KVF_NAMING_H_
