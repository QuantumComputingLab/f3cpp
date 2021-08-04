//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_io_INIFile_hpp
#define f3c_io_INIFile_hpp

#include <memory>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>

namespace f3c {

  namespace io {

    /// INI File Handler
    class INIFile
    {

      public:
        /// Constructs an INI file handler from the given filename.
        INIFile( const std::string filename ) ;

        /// Checks if this INI file handler contains a data field `query`.
        bool contains( const std::string query ) ;

        /// Template function which returns the value of data field `query`.
        template <typename T>
        T value( const std::string query ) ;

        /// Template function which returns the vector of data field `query`.
        template <typename T>
        std::vector< T > vector( const std::string query ) ;

        /// Converts this INI file handler to std::string.
        std::string toString() const ;

      private:
        /// Parses this INI file handler.
        void parse() ;

        /// Splits the given key into name and value.
        int split( const std::string& key ,
                   std::string& name ,
                   std::string& value ,
                   const std::string delimiter ) const ;

        /// Input stream of this INI file handler.
        std::ifstream                                            stream_ ;
        /// Dictionary of this INI file handler.
        std::unordered_map< std::string ,
                            std::unordered_map< std::string ,
                                                std::string > >  dictionary_ ;

    } ; // class INIFile

  } // namespace io

} // namespace f3c

#endif

