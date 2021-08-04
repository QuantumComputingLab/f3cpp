//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_turnoverSU2_hpp
#define f3c_turnoverSU2_hpp

#include "f3c/util.hpp"
#include "qclab/QRotation.hpp"

namespace f3c {

  /// Computes the turnover operation on 3 quantum rotations that form SU(2).
  template <typename T>
  std::tuple< qclab::QRotation< T > ,
              qclab::QRotation< T > ,
              qclab::QRotation< T > >
  turnoverSU2( const qclab::QRotation< T >& rot1 ,
               const qclab::QRotation< T >& rot2 ,
               const qclab::QRotation< T >& rot3 ) {

    const auto rot1times3 = rot1 * rot3 ;
    const auto rot1div3 = rot1 / rot3 ;

    auto ar =  rot2.cos() * rot1times3.cos() ;
    auto ai = -rot2.sin() * rot1div3.cos() ;
    auto br =  rot2.cos() * rot1times3.sin() ;
    auto bi = -rot2.sin() * rot1div3.sin() ;

    // rotation B
    auto cb = std::sqrt( ar*ar + ai*ai ) ;
    auto sb = std::sqrt( br*br + bi*bi ) ;
    qclab::QRotation< T >  rotB( cb , sb ) ;

    if ( cb != 0 ) {
      ar = ar / cb ;
      ai = ai / cb ;
    }
    if ( sb != 0 ) {
      br = br / sb ;
      bi = bi / sb ;
    }

    // rotation A
    T ca ;
    T sa ;
    if ( ( ar + br ) * ( ar + br ) + ( bi - ai ) * ( bi - ai ) >=
         ( ai + bi ) * ( ai + bi ) + ( br - ar ) * ( br - ar ) ) {
      std::tie( ca , sa ) = rotateToZero( ar + br , bi - ai ) ;
    } else {
      std::tie( ca , sa ) = rotateToZero( -ai - bi , br - ar ) ;
    }
    qclab::QRotation< T >  rotA( ca , sa ) ;

    // rotation C
    T cc ;
    T sc ;
    if ( ( bi - ai ) * ( bi - ai ) + ( br - ar ) * ( br - ar )  >
         ( ar + br ) * ( ar + br ) + ( ai + bi ) * ( ai + bi ) ) {
      std::tie( cc , sc ) = rotateToZero( bi - ai , br - ar ) ;
    } else {
      std::tie( cc , sc ) = rotateToZero( ar + br , -ai - bi ) ;
    }
    if ( sign( ar ) != sign( cb * ( ca * cc - sa * sc ) ) ) {
      cc = -cc ;
      sc = -sc ;
    }
    qclab::QRotation< T >  rotC( cc , sc ) ;

    return { rotA , rotB , rotC } ;

  }

} // namespace f3c

#endif

