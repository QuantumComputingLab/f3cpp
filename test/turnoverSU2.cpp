#include <gtest/gtest.h>
#include "f3c/turnoverSU2.hpp"

template <typename T>
void test_f3c_turnoverSU2() {

  const T pi = 4 * std::atan(1) ;
  const T eps = std::numeric_limits< T >::epsilon() ;

  {
    qclab::QRotation< T > rot1 ;
    qclab::QRotation< T > rot2 ;
    qclab::QRotation< T > rot3 ;
    auto [ rotA , rotB , rotC ] = f3c::turnoverSU2( rot1 , rot2 , rot3 ) ;
    EXPECT_NEAR( rotA.theta() , 0 , eps ) ;
    EXPECT_NEAR( rotB.theta() , 0 , eps ) ;
    EXPECT_NEAR( rotC.theta() , 0 , eps ) ;
  }

  // TODO: add more tests

}


/*
 * float
 */
TEST( f3c_turnoverSU2 , float ) {
  test_f3c_turnoverSU2< float >() ;
}

/*
 * double
 */
TEST( f3c_turnoverSU2 , double ) {
  test_f3c_turnoverSU2< double >() ;
}

