//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_turnover_hpp
#define f3c_turnover_hpp

#include "f3c/util.hpp"
#include "f3c/turnoverSU2.hpp"
#include "f3c/qgates/RotationXY.hpp"
#include "f3c/qgates/RotationTFXY.hpp"
#include "f3c/qgates/RotationTFXYMatrix.hpp"

namespace f3c {

  /// Computes the turnover operation of 3 XY-rotation gates.
  template <typename T>
  void turnover( const f3c::qgates::RotationXY< T >& gate1 ,
                 const f3c::qgates::RotationXY< T >& gate2 ,
                 const f3c::qgates::RotationXY< T >& gate3 ,
                 std::unique_ptr< f3c::qgates::RotationXY< T > >& gateA ,
                 std::unique_ptr< f3c::qgates::RotationXY< T > >& gateB ,
                 std::unique_ptr< f3c::qgates::RotationXY< T > >& gateC ) {

    // checks
    const auto q1 = gate1.qubits() ;
    const auto q2 = gate2.qubits() ;
    assert( q1[0] == gate3.qubits()[0] ) ;
    assert( q1[1] == gate3.qubits()[1] ) ;
    assert( ( q2[0] == q1[1] ) || ( q2[1] == q1[0] ) ) ;

    // 2 SU(2) turnovers
    const auto& [ rot10 , rot11 ] = gate1.rotations() ;
    const auto& [ rot20 , rot21 ] = gate2.rotations() ;
    const auto& [ rot30 , rot31 ] = gate3.rotations() ;
    const auto [ rotA0 , rotB0 , rotC0 ] = turnoverSU2( rot10 , rot21 , rot30 );
    const auto [ rotA1 , rotB1 , rotC1 ] = turnoverSU2( rot11 , rot20 , rot31 );

    // new gates
    using RXY = f3c::qgates::RotationXY< T > ;
    gateA = std::make_unique< RXY >( q2[0] , q2[1] , rotA1 , rotA0 ) ;
    gateB = std::make_unique< RXY >( q1[0] , q1[1] , rotB0 , rotB1 ) ;
    gateC = std::make_unique< RXY >( q2[0] , q2[1] , rotC1 , rotC0 ) ;

  }

  /// Computes the turnover operation of 3 TFXY-rotation gates.
  template <typename T>
  void turnover( const f3c::qgates::RotationTFXYMatrix< T >& gate1 ,
                 const f3c::qgates::RotationTFXYMatrix< T >& gate2 ,
                 const f3c::qgates::RotationTFXYMatrix< T >& gate3 ,
              std::unique_ptr< f3c::qgates::RotationTFXYMatrix< T > >& gateA ,
              std::unique_ptr< f3c::qgates::RotationTFXYMatrix< T > >& gateB ,
              std::unique_ptr< f3c::qgates::RotationTFXYMatrix< T > >& gateC ) {

    // checks
    const auto q1 = gate1.qubits() ;
    const auto q2 = gate2.qubits() ;
    assert( q1[0] == gate3.qubits()[0] ) ;
    assert( q1[1] == gate3.qubits()[1] ) ;
    assert( ( q2[0] == q1[1] ) || ( q2[1] == q1[0] ) ) ;

    // work arrays
    std::array< T , 4 >  Q1 ;
    std::array< T , 4 >  Q2 ;
    std::array< T , 4 >  U ;
    std::array< T , 4 >  V ;
    std::array< T , 4 >  Y ;
    std::array< T , 4 >  Z ;
    std::array< T , 4 >  W ;
    T* v ;

    // gate values
    const auto& v1 = gate1.values() ;
    const auto& v2 = gate2.values() ;
    const auto& v3 = gate3.values() ;

    // matrix blocks
    if ( q2[0] > q1[0] ) {
      Q1[0] =  v3[0] * v2[0]            * v1[0] - std::conj(v3[3]) * std::conj(v2[1]) * v1[3] ;
      Q1[1] =  v3[1] * v2[3]            * v1[0] + std::conj(v3[2]) * std::conj(v2[2]) * v1[3] ;
      Q1[2] = -v3[0] * std::conj(v2[3]) * v1[1] - std::conj(v3[3]) * v2[2]            * v1[2] ;
      Q1[3] =  v3[1] * std::conj(v2[0]) * v1[1] - std::conj(v3[2]) * v2[1]            * v1[2] ;

      Q2[0] =  v3[2] * v2[3]            * v1[0] - std::conj(v3[1]) * std::conj(v2[2]) * v1[3] ;
      Q2[1] =  v3[3] * v2[0]            * v1[0] + std::conj(v3[0]) * std::conj(v2[1]) * v1[3] ;
      Q2[2] =  v3[2] * std::conj(v2[0]) * v1[1] + std::conj(v3[1]) * v2[1]            * v1[2] ;
      Q2[3] = -v3[3] * std::conj(v2[3]) * v1[1] + std::conj(v3[0]) * v2[2]            * v1[2] ;
    } else {
      Q1[0] =  v3[0]            * v2[0]            * v1[0]            - std::conj(v3[3]) * v2[1]            * v1[3] ;
      Q1[1] =  std::conj(v3[1]) * v2[3]            * v1[0]            + v3[2]            * v2[2]            * v1[3] ;
      Q1[2] = -v3[0]            * std::conj(v2[3]) * std::conj(v1[1]) - std::conj(v3[3]) * std::conj(v2[2]) * std::conj(v1[2]) ;
      Q1[3] =  std::conj(v3[1]) * std::conj(v2[0]) * std::conj(v1[1]) - v3[2]            * std::conj(v2[1]) * std::conj(v1[2]) ;

      Q2[0] = -std::conj(v3[2]) * v2[3]            * v1[0]            + v3[1]            * v2[2]            * v1[3] ;
      Q2[1] =  v3[3]            * v2[0]            * v1[0]            + std::conj(v3[0]) * v2[1]            * v1[3] ;
      Q2[2] = -std::conj(v3[2]) * std::conj(v2[0]) * std::conj(v1[1]) - v3[1]            * std::conj(v2[1]) * std::conj(v1[2]) ;
      Q2[3] = -v3[3]            * std::conj(v2[3]) * std::conj(v1[1]) + std::conj(v3[0]) * std::conj(v2[2]) * std::conj(v1[2]) ;
    }

    // turnover
    const auto eps = std::numeric_limits< qclab::real_t< T > >::epsilon() ;
    if ( ( norm22( Q1.data() ) - norm22( Q2.data() ) < 100*eps ) ||
         ( std::abs(Q2[1]) + std::abs(Q2[2]) == 0 ) ) {
      //
      // diagonalize, anti-diagonalize, diagonalize, anti-diagonalize
      //

      // (1) compute U, Y that diagonalize (Q11,Q33)
      Z[0] =  std::conj(Q1[3]) ;
      Z[1] = -std::conj(Q1[2]) ;
      Z[2] = -std::conj(Q1[1]) ;
      Z[3] =  std::conj(Q1[0]) ;
      diagonalize22( Q1.data() , Z.data() , U.data() , Y.data() ) ;

      // (2) compute V to anti-diagonalize (Q41,Q23)
      Z[0] =  std::conj(Q2[3]) ;
      Z[1] = -std::conj(Q2[2]) ;
      Z[2] = -std::conj(Q2[1]) ;
      Z[3] =  std::conj(Q2[0]) ;
      gemm22(  Z.data() , Y.data() , W.data() ) ;  // AD23
      gemm22( Q2.data() , Y.data() , Z.data() ) ;  // AD41
      if ( std::abs(Z[1]) + std::abs(Z[3]) > std::abs(W[1]) + std::abs(W[3]) ) {
        // use AD41
        if ( std::abs( Z[3] ) > std::abs( Z[1] ) ) {
          rotateToZeroL( Z[2] , Z[3] , V.data() ) ;
        } else {
          rotateToZeroL( Z[1] , Z[0] , V.data() ) ;
          // V = rot90(V,2)
          std::swap( V[0] , V[3] ) ;
          std::swap( V[1] , V[2] ) ;
        }
      } else {
        // use AD23
        if ( std::abs( W[3] ) > std::abs( W[1] ) ) {
          rotateToZeroL( W[2] , W[3] , V.data() ) ;
        } else {
          rotateToZeroL( W[1] , W[0] , V.data() ) ;
          // V = rot90(V,2)
          std::swap( V[0] , V[3] ) ;
          std::swap( V[1] , V[2] ) ;
        }
      }
      Q1[1] = V[1] * Z[0] + V[3] * Z[1] ;
      Q1[2] = V[0] * Z[2] + V[2] * Z[3] ;

      // (3) compute Z to diagonalize (Q22,Q44)
      if ( q2[0] > q1[0] ) {
        Q2[0] =  v3[0] * v2[1]            * v1[0] - std::conj(v3[3]) * std::conj(v2[0]) * v1[3] ;
        Q2[1] =  v3[1] * v2[2]            * v1[0] + std::conj(v3[2]) * std::conj(v2[3]) * v1[3] ;
        Q2[2] = -v3[0] * std::conj(v2[2]) * v1[1] - std::conj(v3[3]) * v2[3]            * v1[2] ;
        Q2[3] =  v3[1] * std::conj(v2[1]) * v1[1] - std::conj(v3[2]) * v2[0]            * v1[2] ;
      } else {
        Q2[0] =  v3[0]            * std::conj(v2[1]) * v1[0]            - std::conj(v3[3]) * std::conj(v2[0]) * v1[3] ;
        Q2[1] = -std::conj(v3[1]) * std::conj(v2[2]) * v1[0]            - v3[2]            * std::conj(v2[3]) * v1[3] ;
        Q2[2] =  v3[0]            * v2[2]            * std::conj(v1[1]) + std::conj(v3[3]) * v2[3]            * std::conj(v1[2]) ;
        Q2[3] =  std::conj(v3[1]) * v2[1]            * std::conj(v1[1]) - v3[2]            * v2[0]            * std::conj(v1[2]) ;
      }
      rotateToZeroR( V[1] * Q2[2] + V[3] * Q2[3] ,
                     V[1] * Q2[0] + V[3] * Q2[1] , Z.data() ) ;
      rotateToZeroR( -V[1] * std::conj(Q2[1]) + V[3] * std::conj(Q2[0]) ,
                      V[1] * std::conj(Q2[3]) - V[3] * std::conj(Q2[2]) ,
                     W.data() ) ;
      // free Y
      U[1] = Y[0] ;
      U[3] = Y[2] ;
      gemm22( V.data() , Q2.data() , Y.data() ) ;
      if ( ( std::abs( Q1[0] - ( Y[0]*W[0] + Y[2]*W[1] ) ) +
             std::abs( Q1[3] - ( Y[1]*W[2] + Y[3]*W[3] ) ) ) <
           ( std::abs( Q1[0] - ( Y[0]*Z[0] + Y[2]*Z[1] ) ) +
             std::abs( Q1[3] - ( Y[1]*Z[2] + Y[3]*Z[3] ) ) ) ) {
        Z[0] = W[0] ; Z[1] = W[1] ; Z[2] = W[2] ; Z[3] = W[3] ;
      }

      // (4) set pointer to values for gateB
      v = &Q1[0] ;
      Y[0] = U[1] ;
      Y[2] = U[3] ;

    } else {
      //
      // anti-diagonalize, diagonalize, anti-diagonalize, diagonalize
      //

      // (1) compute V, Y that anti-diagonalize (Q23,Q41)
      Z[0] =  std::conj(Q2[3]) ;
      Z[1] =  std::conj(Q2[2]) ;
      Z[2] = -std::conj(Q2[1]) ;
      Z[3] = -std::conj(Q2[0]) ;
      Q2[1] = -Q2[1] ;
      Q2[3] = -Q2[3] ;
      diagonalize22( Z.data() , Q2.data() , V.data() , Y.data() ) ;
      Q2[1] = Q2[0] ;
      Q2[2] = Q2[3] ;
      // V = flipup(V)
      std::swap( V[0] , V[1] ) ;
      std::swap( V[2] , V[3] ) ;
      // V(:,2) = -V(:,2)
      V[2] = -V[2] ;
      V[3] = -V[3] ;

      // (2) compute U to diagonalize (Q11,Q33)
      gemm22( Q1.data() , Y.data() , Z.data() ) ;
      Q2[0] = -std::conj(Q1[2]) * Y[0] + std::conj(Q1[0]) * Y[1] ;
      if ( std::abs( Z[1] ) > std::abs( Q2[0] ) ) {
        rotateToZeroL( Z[0] , Z[1] , U.data() ) ;
      } else {
        rotateToZeroL( std::conj(Q1[3]) * Y[0] - std::conj(Q1[1]) * Y[1] ,
                       Q2[0] , U.data() ) ;
      }
      // D = U * Z
      Q2[0] = U[0] * Z[0] + U[2] * Z[1] ;
      Q2[3] = U[1] * Z[2] + U[3] * Z[3] ;

      // (3) compute Z to anti-diagonalize (Q14,Q32)
      if ( q2[0] > q1[0] ) {
        Q1[0] =  v3[2] * v2[2]            * v1[0] - std::conj(v3[1]) * std::conj(v2[3]) * v1[3] ;
        Q1[1] =  v3[3] * v2[1]            * v1[0] + std::conj(v3[0]) * std::conj(v2[0]) * v1[3] ;
        Q1[2] =  v3[2] * std::conj(v2[1]) * v1[1] + std::conj(v3[1]) * v2[0]            * v1[2] ;
        Q1[3] = -v3[3] * std::conj(v2[2]) * v1[1] + std::conj(v3[0]) * v2[3]            * v1[2] ;
      } else {
        Q1[0] =  std::conj(v3[2]) * std::conj(v2[2]) * v1[0]            - v3[1]            * std::conj(v2[3]) * v1[3] ;
        Q1[1] =  v3[3]            * std::conj(v2[1]) * v1[0]            + std::conj(v3[0]) * std::conj(v2[0]) * v1[3] ;
        Q1[2] = -std::conj(v3[2]) * v2[1]            * std::conj(v1[1]) - v3[1]            * v2[0]            * std::conj(v1[2]) ;
        Q1[3] =  v3[3]            * v2[2]            * std::conj(v1[1]) - std::conj(v3[0]) * v2[3]            * std::conj(v1[2]) ;
      }
      rotateToZeroR( U[0] * Q1[2] + U[2] * Q1[3] ,
                     U[0] * Q1[0] + U[2] * Q1[1] , Z.data() ) ;
      rotateToZeroR( -U[0] * std::conj(Q1[1]) + U[2] * std::conj(Q1[0]) ,
                      U[0] * std::conj(Q1[3]) - U[2] * std::conj(Q1[2]) ,
                     W.data() ) ;
      // free Y
      V[1] = Y[0] ;
      V[3] = Y[2] ;
      gemm22( U.data() , Q1.data() , Y.data() ) ;
      if ( ( std::abs( Q2[1] - ( Y[1]*W[0] + Y[3]*W[1] ) ) +
             std::abs( Q2[2] - ( Y[0]*W[2] + Y[2]*W[3] ) ) ) <
           ( std::abs( Q2[1] - ( Y[1]*Z[0] + Y[3]*Z[1] ) ) +
             std::abs( Q2[2] - ( Y[0]*Z[2] + Y[2]*Z[3] ) ) ) ) {
        Z[0] = W[0] ; Z[1] = W[1] ; Z[2] = W[2] ; Z[3] = W[3] ;
      }

      // (4) set pointer to values for gateB
      v = &Q2[0] ;
      Y[0] = V[1] ;
      Y[2] = V[3] ;

    }

    // new gates
    using TFXY = f3c::qgates::RotationTFXYMatrix< T > ;
    if ( q2[0] > q1[0] ) {
      // vee --> hat
      gateA = std::make_unique< TFXY >( q2[0] , q2[1] , std::conj( Y[0] ) ,
                                                        std::conj( Z[0] ) ,
                                                        std::conj( Z[2] ) ,
                                                        std::conj( Y[2] ) ) ;
      gateB = std::make_unique< TFXY >( q1[0] , q1[1] , v[0] , v[3] ,
                                                        v[2] , v[1] ) ;
      gateC = std::make_unique< TFXY >( q2[0] , q2[1] , std::conj( U[0] ) ,
                                                        std::conj( V[0] ) ,
                                                        std::conj( V[2] ) ,
                                                        std::conj( U[2] ) ) ;
    } else {
      // hat --> vee
      gateA = std::make_unique< TFXY >( q2[0] , q2[1] , std::conj( Y[0] ) ,
                                         Z[0] , -Z[2] , std::conj( Y[2] ) ) ;
      gateB = std::make_unique< TFXY >( q1[0] , q1[1] , v[0] ,
                                                        std::conj( v[3] ) ,
                                                       -std::conj( v[2] ) ,
                                                        v[1] ) ;
      gateC = std::make_unique< TFXY >( q2[0] , q2[1] , std::conj( U[0] ) ,
                                         V[0] , -V[2] , std::conj( U[2] ) ) ;
    }

  }

  /// Computes the turnover operation of 3 gates.
  template <typename G>
  std::tuple< G , G , G > turnover( const G& gate1 ,
                                    const G& gate2 ,
                                    const G& gate3 ) {
    std::unique_ptr< G >  gateA ;
    std::unique_ptr< G >  gateB ;
    std::unique_ptr< G >  gateC ;
    turnover( gate1 , gate2 , gate3 , gateA , gateB , gateC ) ;
    return { *gateA , *gateB , *gateC } ;
  }

} // namespace f3c

#endif

