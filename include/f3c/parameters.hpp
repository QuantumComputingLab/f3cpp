//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_parameters_hpp
#define f3c_parameters_hpp

#include <cassert>
#include <limits>
#include <vector>
#include <sstream>

namespace f3c {

  template <typename R>
  class Param
  {

    public:
      using value_type = R ;

      virtual R value( const size_t timestep ) const = 0 ;

      virtual size_t size() const = 0 ;

      const R operator[]( size_t timestep ) const {
        return value( timestep ) ;
      }

      virtual std::string toString() const = 0 ;

      virtual ~Param() noexcept = default ;

  } ;

  template <typename R>
  std::ostream& operator<<( std::ostream& os , const Param< R >& obj ) {
    return os << obj.toString() ;
  }


  template <typename R>
  class ConstValue : public Param< R >
  {

    public:
      using value_type = R ;

      ConstValue()
      : value_( 0 )
      { }

      ConstValue( const R value )
      : value_( value )
      { }

      R value() const { return value_ ; }

      R value( const size_t timestep ) const override {
        return value_ ;
      }

      size_t size() const override { return 1 ; }

      virtual std::string toString() const override {
        std::stringstream stream ;
        stream << "const value: " << value_ ;
        return stream.str() ;
      }

    protected:
      R value_ ;

  } ;


  template <typename R>
  class Values : public Param< R >
  {

    public:
      using value_type = R ;

      Values()
      { }

      Values( const R value )
      : values_( { value } )
      { }

      Values( const std::vector< R >& values )
      : values_( values )
      { }

      std::vector< R >& values() const { return values_ ; }

      R value( const size_t timestep ) const override {
        assert( timestep < values_.size() ) ;
        return values_[ timestep ] ;
      }

      size_t size() const override { return values_.size() ; }

      virtual std::string toString() const override {
        std::stringstream stream ;
        for ( size_t i = 0; i < values_.size(); i++ ) {
          if ( i > 0 ) stream << "," ;
          stream << values_[i] ;
        }
        return stream.str() ;
      }

    protected:
      std::vector< R >  values_ ;

  } ;


  template <typename R>
  class LinearRamp : public Param< R >
  {

    public:
      using value_type = R ;

      LinearRamp()
      : value_( 0 )
      , level_( 0 )
      , begin_( 0 )
      , end_( std::numeric_limits< size_t >::max() )
      { }

      LinearRamp( const R beginLevel , const R endLevel )
      : value_( beginLevel )
      , level_( endLevel )
      , begin_( 0 )
      , end_( std::numeric_limits< size_t >::max() )
      { }

      LinearRamp( const R beginLevel , const R endLevel ,
                  const size_t beginRampStep , const size_t endRampStep )
      : value_( beginLevel )
      , level_( endLevel )
      , begin_( beginRampStep )
      , end_( endRampStep )
      {
        assert( beginRampStep <= endRampStep ) ;
      }

      R beginLevel()    const { return value_ ; } ;

      R endLevel()      const { return level_ ; } ;

      R beginRampStep() const { return begin_ ; } ;

      R endRampStep()   const { return end_ ; } ;

      R value( const size_t timestep ) const override {
        if ( timestep <= begin_ ) {
          return value_ ;
        } else if ( timestep >= end_ ) {
          return level_ ;
        } else {
          const auto diff = level_ - value_ ;
          return value_ + R(timestep - begin_) / R(end_ - begin_) * diff ;
        }
      }

      size_t size() const override { return 0 ; }

      virtual std::string toString() const override {
        std::stringstream stream ;
        stream << "linear ramp: from " << value_ << " to " << level_ <<
                                  " [" << begin_+1 << ":" << end_+1 << "]" ;
        return stream.str() ;
      }

    protected:
      R       value_ ;
      R       level_ ;
      size_t  begin_ ;
      size_t  end_ ;

  } ;

} // namespace f3c

#endif

