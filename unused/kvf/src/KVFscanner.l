%{

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>

#include "KVFparser.h"

namespace KVF {
#include "KVFyacc.h"
}

long ascii_to_long( long base, const char *s )
{
  long l = 0;
  while( *s != '\0' ) {
    char c = *s++;
    if( c >= 'a' )
      c -= 'a' - 'A';
    c -= '0';
    if( c > 9 )
      c -= 'A' - '0' - 10;
    l = l * base + c;
  }
  return l;
}


char translate_char( const char *s )
{
  char c = *s++;
  
  if( c != '\\' )
    return c;
  c = *s++;
  switch( c ) {
  case 'n':
    return '\n';
  case 't':
    return '\t';
  case 'v':
    return '\v';
  case 'b':
    return '\b';
  case 'r':
    return '\r';
  case 'f':
    return '\f';
  case 'a':
    return '\a';
  case '\\':
    return '\\';
  case '?':
    return '\?';
  case '\'':
    return '\'';
  case '"':
    return '"';
  case 'x':
  case 'X':
    return (char) ascii_to_long( 16, s );
  default:
    // Gotta be an octal
    return (char) ascii_to_long( 8, s );
  }
}


%}

%option c++

%option noyywrap

/*--------------------------------------------------------------------------*/

Digits                  [0-9]+
Oct_Digit               [0-7]
Hex_Digit               [a-fA-F0-9]
Int_Literal		("+"|"-")?[1-9][0-9]*
Oct_Literal		0{Oct_Digit}*
Hex_Literal		(0x|0X){Hex_Digit}*
Esc_Sequence1           "\\"[ntvbrfa\\\?\'\"]
Esc_Sequence2           "\\"{Oct_Digit}{1,3}
Esc_Sequence3           "\\"(x|X){Hex_Digit}{1,2}
Esc_Sequence            ({Esc_Sequence1}|{Esc_Sequence2}|{Esc_Sequence3})
Char                    ([^\n\t\"\'\\]|{Esc_Sequence})
Char_Literal            "'"({Char}|\")"'"
String_Literal		\"({Char}|"'")*\"
Float_Literal1		("+"|"-")?{Digits}"."{Digits}(e|E)("+"|"-")?{Digits}  
Float_Literal2		("+"|"-")?{Digits}(e|E)("+"|"-")?{Digits}
Float_Literal3          ("+"|"-")?{Digits}"."{Digits}
Float_Literal4		("+"|"-")?"."{Digits} 
Float_Literal5		("+"|"-")?"."{Digits}(e|E)("+"|"-")?{Digits}  

/*--------------------------------------------------------------------------*/

KVF_Identifier	[a-zA-Z/][a-zA-Z0-9_/]*

/*--------------------------------------------------------------------------*/



%%

[ \t]			;
[\n]			KVF_line_no++;
"//"[^\n]*		;
"#"[^\n]*		;
 
"{"			return T_LEFT_CURLY_BRACKET;
"}"			return T_RIGHT_CURLY_BRACKET;
"["		 	return T_LEFT_SQUARE_BRACKET;
"]"			return T_RIGHT_SQUARE_BRACKET;
";"			return T_SEMI_COLON;
"="			return T_EQUALS;
","			return T_COMMA;

{KVF_Identifier}	{
			  KVFyylval.ident = strdup( YYText() );
			  return T_IDENTIFIER;
			}
{Float_Literal1}	|
{Float_Literal2}	|
{Float_Literal3}	|
{Float_Literal4}	|
{Float_Literal5}	{
			  KVFyylval._float = atof( YYText() );
			  return T_FLOATING_PT_LITERAL;
			}
{Int_Literal}		{ 
			 /* printf("Int_Literal: %s\n", YYText() ); */
			  KVFyylval._int = strtol( YYText(), NULL, 0 );
			  return T_INTEGER_LITERAL;
			}
{Oct_Literal}		{
			  KVFyylval._int = strtol( YYText(), NULL, 0 );
			  return T_INTEGER_LITERAL;
			}
{Hex_Literal}		{
			  KVFyylval._int = strtol( YYText(), NULL, 0 );
			  return T_INTEGER_LITERAL;
			}
{Char_Literal}		{
                          string s( YYText() );
                          s = s.substr( 1, s.length() - 2 );
			  KVFyylval._char = translate_char( s.c_str() );
			  return T_CHARACTER_LITERAL;
			}
{String_Literal}	{
                          string s( YYText() );
                          s = s.substr( 1, s.length() - 2 );
			  KVFyylval._str = strdup( s.c_str() );
			  return T_STRING_LITERAL;
			}


%%


