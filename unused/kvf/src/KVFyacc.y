
%{
/*#include <iostream> */
#include <stdio.h>
#include <stdlib.h>
#include "KVFparser.h"

namespace KVF {

void KVFyyerror(const char *s)
{
// printf("KVFyyerror: line: %i %s \n", KVF_line_no, s );
 theParser->set_parse_error( s, KVF_line_no );
}

%}

%union 
{
  char         *ident;
  int            _int;
  char         *_str;
  char           _char;
  double         _float;
}


%token <_char> T_CHARACTER_LITERAL
%token <_int> T_INTEGER_LITERAL
%token <_float> T_FLOATING_PT_LITERAL
%token <_str> T_STRING_LITERAL

%token <ident> T_IDENTIFIER

%token T_LEFT_CURLY_BRACKET
%token T_LEFT_PARANTHESIS
%token T_LEFT_SQUARE_BRACKET

%token T_RIGHT_CURLY_BRACKET
%token T_RIGHT_PARANTHESIS
%token T_RIGHT_SQUARE_BRACKET
%token T_SEMI_COLON
%token T_EQUALS
%token T_COMMA

%type <ident> key
%type <_str> string
%type <_float> number

%%

kv_items: /* empty */ {;}
        | kv_item kv_items
        | kv_item
        ;


key:     T_IDENTIFIER { $$ = $1; }
        ;

kv_item:  kv_string
        | kv_str_seq
        | kv_numeric
        | kv_numeric_seq
        ;


kv_string: key T_EQUALS string  T_SEMI_COLON
        { /* printf("kv_string %s %s\n", $1, $3 ); */
          theParser->assign_string( $1, $3 ); }
        ;

kv_numeric: key T_EQUALS number T_SEMI_COLON
        { /* printf("kv_numeric(number) %s %g\n", $1, $3 ); */
          theParser->assign_number( $1, $3 );}
        ;

string: T_STRING_LITERAL { $$=$1; }
        ;

string_seq: string T_COMMA string_seq
          { /* printf("string_seq: seq: %s\n",$1); */
            theParser->push_strseq($1);
          }
        | string
          {
            /* printf("string_seq: single: %s\n",$1); */
            theParser->start_strseq($1);
          }
        ;

kv_str_seq: key T_EQUALS T_LEFT_CURLY_BRACKET string_seq
                            T_RIGHT_CURLY_BRACKET T_SEMI_COLON
           { /* printf("kv_string_seq: %s\n", $1 );
             theParser->diag_print_strseq();  */
             theParser->assign_strseq( $1 ); }
        ;

number:   T_FLOATING_PT_LITERAL { $$=$1; }
        | T_INTEGER_LITERAL { $$=$1; }
        ;

number_seq: number T_COMMA number_seq
          { /* printf("number_seq: seq: %g\n",$1); */
            theParser->push_nseq($1);
          }
        | number
          { /* printf("number_seq: single: %g\n",$1); */
            theParser->start_nseq($1);
          }
        ;       

kv_numeric_seq: key T_EQUALS T_LEFT_CURLY_BRACKET number_seq
                             T_RIGHT_CURLY_BRACKET T_SEMI_COLON
           { /* printf("kv_numeric_seq: %s\n", $1 );
             theParser->diag_print_numseq();  */
             theParser->assign_nseq( $1 ); }
        ;

%%
} // end namespace KVF
