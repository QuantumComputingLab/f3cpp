//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_io_param_hpp
#define f3c_io_param_hpp

#include "f3c/parameters.hpp"
#include "f3c/io/INIFile.hpp"

namespace f3c::io {

  template <typename R>
  int param( const int N , const int n , INIFile& file ,
             const std::string name , std::unique_ptr< Param< R > >& ptr ) {
    if ( file.contains( name + ".value" ) ) {
      ptr = std::make_unique< ConstValue< R > >(
              file.value< double >( name + ".value" ) ) ;
    } else if ( file.contains( name + ".values" ) ) {
      ptr = std::make_unique< Values< R > >(
              file.vector< double >( name + ".values" ) ) ;
      if ( ptr->size() != n ) {
        std::cout << "ERROR: number of values for parameter " << name
                  << " should be equal to " << n << "!" << std::endl ;
        return -2 ;
      }
    } else if ( file.contains( name + ".ramp" ) ) {
      ptr = std::make_unique< LinearRamp< R > >(
              file.value< double >( name + ".value1" ) ,
              file.value< double >( name + ".value2" ) ,
              file.value< int >(    name + ".begin" ) - 1 ,
              file.value< int >(    name + ".end" )   - 1 ) ;
    } else {
      std::cout << "ERROR: parameter " << name << " is required!" << std::endl ;
      return -1 ;
    }
    return 0 ;
  }

} // namespace f3c::io

#endif

