// 
// KVFDataSource.h
// 
// Key Value File data source
// 
// David Burgess
// March 2000, June 2004, September 2006
// 

#ifndef KVF_KVFDATASOURCE_H
#define KVF_KVFDATASOURCE_H


#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "kvfException.h"
#include "kvfNaming.h"
#include "kvfDataBuffer.h"

// extern int KVFyyparse(); // Defined through yacc.y

namespace KVF {

  using namespace std;

  typedef NamingContext < string > DataNamesTable;
  typedef Var < DataNamesTable > DataNamesTable_var;

  /* ! \brief DataSource with data from a KVF (Key-Value-Format) file.
   * 
   * A KVFDataSource opens and reads a KVF file, and stores all of the data in the file, arranged by the key name for each data item.
   * 
   * 
   * \note Current strategy is open, read, and close the file in one operation, either at object construction, or by using read_file().
   * 
   * \todo Add refresh() to re-read file if it has changed
   * 
   */

  class KVFDataSource {
    string _file_name;
    ifstream _istr;
      NamingContext < DataBuffer_var > _data_by_key;
    DataNamesTable_var _data_names;
      vector < double >_nseq_buff;
      vector < string > _strseq_buff;
    typedef enum KVFparse_state_e { PARSE_NIL, PARSE_NSEQ, PARSE_STRSEQ } KVFParseState;
    KVFParseState _parse_state;
      vector < string > _parse_errors;

    void parse(void);
    void assign_number(const char *key, double x);
    void assign_string(const char *key, const char *s);
    void start_nseq(double x);
    void push_nseq(double x);
    void assign_nseq(const char *key);
    void start_strseq(const char *s);
    void push_strseq(const char *s);
    void assign_strseq(const char *key);
    void set_parse_error(string s, int line_no);
    bool key_is_ok(string key);

    friend int KVFyyparse();
    friend void KVFyyerror(const char *s);

  public:
      KVFDataSource(void) {;
    } KVFDataSource(const string & file_name);
     ~KVFDataSource(void);

    void read_file(const string & file_name);
    void open(const string & file_name);
    void close();

    DataBuffer_var get_data(const string & data_name);

    DataNamesTable_var get_data_names(void);
    void get_data_names(vector < string > &names_list, const string & prefix = "", const string & separator = "/");

    bool name_exists(const string & data_name);

    void get_data(const string & data_name, double &d_val);
    void get_data(const string & data_name, string & s_val);
    void get_data(const string & data_name, vector < double >&d_vvals);
    void get_data(const string & data_name, vector < string > &s_vvals);

    void get_data(const string & data_name, int &i_val);
    void get_data(const string & data_name, vector < int >&i_vvals);
    void get_data(const string & data_name, float &f_val);
    void get_data(const string & data_name, vector < float >&f_vvals);

    int num_parse_errors(void) {
      return _parse_errors.size();
    } void diag_print_strseq() {
      for (int i = 0; i < _strseq_buff.size(); i++)
        cout << _strseq_buff[i] << " ";
      cout << endl;
    }
    void diag_print_numseq() {
      for (int i = 0; i < _nseq_buff.size(); i++)
        cout << _nseq_buff[i] << " ";
      cout << endl;
    }
    void diag_print_errors() {
      for (int i = 0; i < _parse_errors.size(); i++)
        cout << _parse_errors[i] << endl;
    }

    void diag_print_data(void);

  };                            // end class KVFDataSource

}                               // end namespace KVF

#endif // #ifndef KVF_KVFDATASOURCE_H
