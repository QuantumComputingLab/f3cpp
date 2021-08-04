//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_functors_hpp
#define f3c_qgates_functors_hpp

#include "f3c/qgates/RotationXY.hpp"
#include "f3c/qgates/RotationTFXY.hpp"
#include "f3c/qgates/RotationTFXYMatrix.hpp"
#include <memory>

namespace f3c {

  namespace qgates {

    /// XY functor.
    template <typename T>
    struct XYfunctor {

      /// Value type of this XY functor.
      using value_type = T ;
      /// Gate type of this XY functor.
      using gate_type = f3c::qgates::RotationXY< T > ;

      /// Returns a unique pointer to a random XY gate.
      template <typename D, typename G>
      constexpr static inline
      std::unique_ptr< gate_type > init( const int q , D& dis , G& gen ) {
        return std::make_unique< gate_type >( q , q+1 , dis(gen) , dis(gen) ) ;
      }

      /// Constructs 1 timestep with the given parameters.
      template <typename R, typename C>
      constexpr static
      void timestep( const R dt , const R hx , const R hy , const R hz ,
                     const R Jx , const R Jy , const R Jz , C& circuit ) {
        assert( dt > 0 ) ;  assert( Jz == 0 ) ;
        assert( hx == 0 ) ; assert( hy == 0 ) ; assert( hz == 0 ) ;
        const int n = circuit.nbQubits() ;
        assert( circuit.nbGates() == n - 1 ) ;
        // angles
        const auto tJx = 2*dt*Jx ;
        const auto tJy = 2*dt*Jy ;
        // 1st layer
        #pragma omp parallel for
        for ( int i = 0; i < n/2; i++ ) {
          const int q = 2*i ;
          circuit[i] = std::make_unique< gate_type >( q , q+1 , tJx , tJy ) ;
        }
        // 2nd layer
        #pragma omp parallel for
        for ( int i = 0; i < (n-1)/2; i++ ) {
          const int q = 2*i + 1 ;
          circuit[n/2+i] = std::make_unique< gate_type >( q , q+1 , tJx , tJy );
        }
      }

    } ; // XYfunctor


    /// TFXY functor.
    template <typename T>
    struct TFXYfunctor {

      /// Value type of this TFXY functor.
      using value_type = T ;
      /// Gate type of this TFXY functor.
      using gate_type = f3c::qgates::RotationTFXYMatrix< T > ;

      /// Returns a unique pointer to a random TFXY gate.
      template <typename D, typename G>
      constexpr static inline
      std::unique_ptr< gate_type > init( const int q , D& dis , G& gen ) {
        auto r1 = dis( gen ) ;
        auto r2 = std::sqrt( 1 - r1*r1 ) ;
        auto theta1 = dis( gen ) ;
        auto theta2 = dis( gen ) ;
        const T a( r1 * std::cos( theta1 ) , r1 * std::sin( theta1 ) ) ;
        const T d( r2 * std::cos( theta2 ) , r2 * std::sin( theta2 ) ) ;
        r1 = dis( gen ) ;
        r2 = std::sqrt( 1 - r1*r1 ) ;
        theta1 = dis( gen ) ;
        theta2 = dis( gen ) ;
        const T b( r1 * std::cos( theta1 ) , r1 * std::sin( theta1 ) ) ;
        const T c( r2 * std::cos( theta2 ) , r2 * std::sin( theta2 ) ) ;
        return std::make_unique< gate_type >( q , q+1 , a , b , c , d ) ;
      }

      /// Constructs 1 timestep with the given parameters.
      template <typename R, typename C>
      constexpr static
      void timestep( const R dt , const R hx , const R hy , const R hz ,
                     const R Jx , const R Jy , const R Jz , C& circuit ) {
        assert( dt > 0 ) ;  assert( Jz == 0 ) ;
        assert( hx == 0 ) ; assert( hy == 0 ) ;
        const int n = circuit.nbQubits() ;
        assert( circuit.nbGates() == n - 1 ) ;
        // angles
        const auto tJx = 2*dt*Jx ;
        const auto tJy = 2*dt*Jy ;
        const auto thz = 2*dt*hz ;
        // 1st layer
        #pragma omp parallel for
        for ( int i = 0; i < n/2; i++ ) {
          const int q = 2*i ;
          circuit[i] = std::make_unique< gate_type >( q , q+1 ,
                                               thz , thz , tJx , tJy , 0 , 0 ) ;
        }
        // 2nd layer
        #pragma omp parallel for
        for ( int i = 0; i < n/2-1; i++ ) {
          const int q = 2*i + 1 ;
          circuit[n/2+i] = std::make_unique< gate_type >( q , q+1 ,
                                               0   , 0   , tJx , tJy , 0 , 0 ) ;
        }
        if ( n % 2 == 1 ) {
          const int q = n - 2 ;
          circuit[n - 2] = std::make_unique< gate_type >( q , q+1 ,
                                               0   , thz , tJx , tJy , 0 , 0 ) ;
        }
      }

    } ; // TFXYfunctor

  } // namespace qgates

} // namespace f3c

#endif

