
#ifndef __STRINGCHOP_H__
#define __STRINGCHOP_H__

namespace std_plus {

     template <class T>
     string vector_join( const vector<T>& v, const string& token ) {
          ostringstream result;
          for (typename vector<T>::const_iterator i = v.begin(); i != v.end(); i++) {
               if (i != v.begin()) { result << token; }
               result << *i;
          }
          return result.str();
     }

     vector<string> string_split( const string& s, const string& delimiter ) {
          vector<string> result;
          string::size_type from = 0;
          string::size_type to = 0;

          while ( to != string::npos ) {
               to = s.find( delimiter, from );
               if ( from < s.size() && from != to ) {
                    result.push_back( s.substr( from, to - from ) );
               }
               from = to + delimiter.size();
          }
          return result;
     }
}

#endif