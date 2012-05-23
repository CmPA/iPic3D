//
//  Naming.h
//
//  Qdos2 - Name Context
//
//  David Burgess
//  June 1999, June 2004, September 2006
//

#include <string>
#include <vector>

#include "kvfNaming.h"

namespace KVF {

using namespace std;

NameComponent& Name::get_component( int i )
{
  if( i < 0 || i >= components.size() )
    throw KVFNamingError("No name components in name", 
            "Name::get_component()" );
  return components[i];
}

const NameComponent& Name::get_component( int i ) const
{
  if( i < 0 || i >= components.size() )
    throw KVFNamingError("No name components in name", 
            "Name::get_component()const" );
  return components[i];
}

NameComponent& Name::operator[]( int i )
{
  return get_component(i);
}

const NameComponent& Name::operator[]( int i ) const
{
  return get_component(i);
}


Name::Name( const string& name, string separators )
{
  int i1=0, i2;
  while( (i2=name.find_first_of( separators, i1 )) !=string::npos ) {
    components.push_back( NameComponent(name.substr(i1,i2-i1)) );
    i1=i2+1;
  }
  if( i1 < name.size() )
    components.push_back( NameComponent(name.substr(i1,name.size()-i1)) );
}

Name::Name( const Name& name )
{
  for(int i=0; i<name.components.size(); i++ )
    components.push_back( name.components[i] );
}

Name::Name( const NameComponent& name_component )
{
  components.push_back( name_component );
}

string Name::srep( const string& separator ) const
{
  string s;
  if( components.size() == 0)
    return s;
  s = components[0].identifier;
  for( int i=1; i<components.size(); ++i )
    s += separator + components[i].identifier;
 return s;
}

void Name::delete_component( int i )
{
  if( i < 0 || i >= components.size() )
    throw KVFNamingError("No name components in name", 
            "Name::delete_component()" );
  components.erase( components.begin() + i );
}

} // end namespace KVF


