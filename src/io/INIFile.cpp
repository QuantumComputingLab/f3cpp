#include "f3c/io/INIFile.hpp"
#include <cassert>
#include <sstream>

namespace f3c::io {

  /// Trim a string of left trailing whitespaces.
  static inline std::string& trim_left( std::string& str ) {
    const auto pos = str.find_first_not_of( " \t" ) ;
    if ( pos == std::string::npos ) {
      str = "" ;
    } else if ( pos > 0 ) {
      str = str.substr( pos ) ;
    }
    return str ;
  } ; // trim_left

  /// Trim a string of right trailing whitespaces.
  static inline std::string& trim_right( std::string& str ) {
    const auto pos = str.find_last_not_of( " \t" ) ;
    if ( pos == std::string::npos ) {
      str = "" ;
    } else if ( pos < str.length() ) {
      str = str.substr( 0 , pos + 1 ) ;
    }
    return str ;
  } ; // trim_right

  /// Trim a string of left and right trailing whitespaces.
  static inline std::string& trim( std::string& str ) {
    return trim_left( trim_right( str ) ) ;
  } ; // trim


  /// Constructs an INI file handler from the given filename.
  INIFile::INIFile( const std::string filename )
  : stream_( std::ifstream( filename ) )
  {
    assert( stream_.good() ) ;
    parse() ;
  } // INIFile(filename)


  /// Checks if this INI file handler contains a data field `query`.
  bool INIFile::contains( const std::string query ) {
    std::string section ;
    std::string field ;
    if ( split( query , section , field , "." ) != 0 ) return false ;
    if ( dictionary_.find( section ) == dictionary_.end() ) return false ;
    return ( dictionary_[ section ].find( field ) !=
             dictionary_[ section ].end() ) ;
  }


  /// Returns std::string of query data field.
  template <>
  std::string INIFile::value( const std::string query ) {
    std::string section ;
    std::string field ;
    if ( split( query , section , field , "." ) == 0 ) {
      return dictionary_[ section ][ field ] ;
    }
    return "" ;
  }

  /// Returns int of query data field.
  template <>
  int INIFile::value( const std::string query ) {
    return std::stoi( value< std::string >( query ) ) ;
  }

  /// Returns float of query data field.
  template <>
  float INIFile::value( std::string query ) {
    return std::stof( value< std::string >( query ) ) ;
  }

  /// Returns double of query data field.
  template <>
  double INIFile::value( std::string query ) {
    return std::stod( value< std::string >( query ) ) ;
  }


  /// Returns std::vector<int> of query data field.
  template <>
  std::vector< std::string > INIFile::vector( std::string query ) {

    // data field string
    std::string str = value< std::string >( query ) ;

    // positions of "[" and "]"
    auto left  = str.find( '[' ) ;
    auto right = str.find( ']' ) ;

    // vector?
    bool vectorField = ( left == 0 ) && ( right == str.length() - 1 ) ;

    // parse
    std::vector< std::string > vec ;
    if ( vectorField ) {
      // vector field
      str = str.substr( 1 , str.length() - 2 ) ;
      std::string val1 ;
      std::string val2 ;
      size_t c = 0 ;
      while ( split( str , val1 , val2 , "," ) == 0 ) {
        vec.push_back( val1 ) ;
        str = val2 ;
        c++ ;
      }
      if ( c == 0 ) {
        // only 1 element
        if (  str.length() > 0 ) { vec.push_back( str ) ; }
      } else {
        // add last element
        if ( val2.length() > 0 ) { vec.push_back( val2 ) ; }
      }
    } else {
      // value field
      vec.push_back( str ) ;
    }
    return vec ;
  }

  /// returns std::vector<int> of query data field.
  template <>
  std::vector< int > INIFile::vector( std::string query ) {
    auto strvec = vector< std::string >( query ) ;
    std::vector< int > vec( strvec.size() ) ;
    #pragma omp parallel for
    for ( size_t i = 0; i < strvec.size(); i++ ) {
      vec[i] = std::stoi( strvec[i] ) ;
    }
    return vec ;
  }

  /// returns std::vector<float> of query data field.
  template <>
  std::vector< float > INIFile::vector( std::string query ) {
    auto strvec = vector< std::string >( query ) ;
    std::vector< float > vec( strvec.size() ) ;
    #pragma omp parallel for
    for ( size_t i = 0; i < strvec.size(); i++ ) {
      vec[i] = std::stof( strvec[i] ) ;
    }
    return vec ;
  }

  /// returns std::vector<double> of query data field.
  template <>
  std::vector< double > INIFile::vector( std::string query ) {
    auto strvec = vector< std::string >( query ) ;
    std::vector< double > vec( strvec.size() ) ;
    #pragma omp parallel for
    for ( size_t i = 0; i < strvec.size(); i++ ) {
      vec[i] = std::stod( strvec[i] ) ;
    }
    return vec ;
  }


  /// Converts this INI file handler to std::string.
  std::string INIFile::toString() const {
    std::stringstream stream ;
    for ( const auto& s : dictionary_ ) {
      stream << "[" << s.first << "]" << std::endl ;
      for ( const auto& k : s.second ) {
        stream << "  " << k.first << "=" << k.second << std::endl ;
      }
    }
    return stream.str() ;
  } // toString()


  /// Parses this INI file handler.
  void INIFile::parse() {

    bool parseSection = false ;
    std::string section ;
    std::string name ;
    std::string value ;

    // loop over lines in file
    while ( not stream_.eof() ) {

      std::string line ;
      std::getline( stream_ , line ) ;

      // skip blank lines
      if ( line.length() < 1 ) { continue ; }

      // positions of first and last non-space character
      auto firstNonSpace = line.find_first_not_of( " " ) ;
      auto lastNonSpace  = line.find_last_not_of( " " ) ;

      // skip comment (lines)
      auto comment = line.find( "#" ) ;
      if ( firstNonSpace == comment ) { continue ; }
      line = line.substr( 0 , comment ) ;

      // strip leading and trailing spaces
      trim( line ) ;

      // positions of "[" , "]" , "=", and ":"
      auto left  = line.find( '[' ) ;
      auto right = line.find( ']' ) ;
      auto equal = line.find( '=' ) ;
      auto colon = line.find( ':' ) ;

      // section header?
      bool sectionLine = ( left == firstNonSpace ) &&
                         ( right == lastNonSpace ) ;

      // data field?
      bool dataLine = ( equal != std::string::npos ) ||
                      ( colon != std::string::npos ) ;

      // section line
      if ( sectionLine ) {
        section = line.substr( 1 , line.length() - 2 ) ;
        dictionary_[ section ] = std::unordered_map< std::string ,
                                                     std::string >() ;
        parseSection = true ;
        continue ;
      }

      // data line
      if ( parseSection && dataLine ) {
        if ( split( line , name , value , "=:" ) == 0 ) {
          dictionary_[ section ][ name ] = value ;
        }
      }

    }

  } // parse()


  /// Splits the given key into name and value.
  int INIFile::split( const std::string& key , std::string& name ,
                      std::string& value , const std::string delimiter ) const {
    auto pos = key.find_first_of( delimiter ) ;
    if ( ( pos == 0 ) || ( pos == std::string::npos ) ) {
      return -1 ;
    } else {
      name = key.substr( 0 , pos ) ;  trim( name ) ;
      value = key.substr( pos + 1 ) ; trim( value ) ;
    }
    return 0 ;
  } // split

}

