#include <gtest/gtest.h>
#include "f3c/SquareCircuit.hpp"
#include "f3c/TriangleCircuit.hpp"
#include "qclab/qgates/RotationXX.hpp"
#include "f3c/qgates/RotationXY.hpp"

template <typename T>
auto test_f3c_SquareCircuit_init( const int n ) {

  using XX = qclab::qgates::RotationXX< T > ;

  f3c::SquareCircuit< T , XX >  square( n ) ;
  int c = 0 ;
  for ( int l = 0; l < n; l++ ) {
    const int qstart = ( l % 2 == 0 ) ? 0 : 1 ;
    for ( int q = qstart; q < n-1; q += 2 ) {
      square[c] = std::make_unique< XX >( q , q+1 , 1./(c+1) ) ;
      c++ ;
    }
  }
  return square ;

}


template <typename T>
auto test_f3c_SquareCircuit_initXY( const int n ) {

  using XY = f3c::qgates::RotationXY< T > ;

  f3c::SquareCircuit< T , XY >  square( n ) ;
  int c = 0 ;
  for ( int l = 0; l < n; l++ ) {
    const int qstart = ( l % 2 == 0 ) ? 0 : 1 ;
    for ( int q = qstart; q < n-1; q += 2 ) {
      square[c] = std::make_unique< XY >( q , q+1 , 1./(c+1) , 1. + 1./(c+1) ) ;
      c++ ;
    }
  }
  return square ;

}


template <typename T>
void test_f3c_SquareCircuit() {

  using R = qclab::real_t< T > ;
  using XY = f3c::qgates::RotationXY< T > ;

  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    f3c::SquareCircuit< T >  square2( 2 ) ;
    EXPECT_EQ( square2.nbQubits() , 2 ) ;
    EXPECT_EQ( square2.capacity() , 1 ) ;

    f3c::SquareCircuit< T >  square3( 3 ) ;
    EXPECT_EQ( square3.nbQubits() , 3 ) ;
    EXPECT_EQ( square3.capacity() , 3 ) ;

    f3c::SquareCircuit< T >  square4( 4 ) ;
    EXPECT_EQ( square4.nbQubits() , 4 ) ;
    EXPECT_EQ( square4.capacity() , 6 ) ;

    f3c::SquareCircuit< T >  square5( 5 ) ;
    EXPECT_EQ( square5.nbQubits() , 5 ) ;
    EXPECT_EQ( square5.capacity() , 10 ) ;

    f3c::SquareCircuit< T >  square6( 6 ) ;
    EXPECT_EQ( square6.nbQubits() , 6 ) ;
    EXPECT_EQ( square6.capacity() , 15 ) ;

    f3c::SquareCircuit< T >  square7( 7 ) ;
    EXPECT_EQ( square7.nbQubits() , 7 ) ;
    EXPECT_EQ( square7.capacity() , 21 ) ;

    f3c::SquareCircuit< T >  square8( 8 ) ;
    EXPECT_EQ( square8.nbQubits() , 8 ) ;
    EXPECT_EQ( square8.capacity() , 28 ) ;
  }

  {
    f3c::SquareCircuit< T , XY >  square4( 4 ) ;
    int c = 1 ;
    for ( auto it = square4.begin(); it != square4.end(); ++it ) {
      *it = std::make_unique< XY >( 1./c , 1. + 1./c ) ;
      c++ ;
    }

    for ( int i = 0; i < square4.nbGates(); i++ ) {
      const int c = i + 1 ;
      const auto [ theta0 , theta1 ] = square4[i]->thetas() ;
      EXPECT_NEAR( theta0 , 1./c , 10*eps ) ;
      EXPECT_NEAR( theta1 , 1. + 1./c , 10*eps ) ;
    }

    auto square5 = test_f3c_SquareCircuit_init< T >( 5 ) ;

    EXPECT_EQ( square5[0]->qubit() , 0 ) ;  // 01
    EXPECT_EQ( square5[1]->qubit() , 2 ) ;  // 02

    EXPECT_EQ( square5[2]->qubit() , 1 ) ;  // 03
    EXPECT_EQ( square5[3]->qubit() , 3 ) ;  // 04

    EXPECT_EQ( square5[4]->qubit() , 0 ) ;  // 05
    EXPECT_EQ( square5[5]->qubit() , 2 ) ;  // 06

    EXPECT_EQ( square5[6]->qubit() , 1 ) ;  // 07
    EXPECT_EQ( square5[7]->qubit() , 3 ) ;  // 08

    EXPECT_EQ( square5[8]->qubit() , 0 ) ;  // 09
    EXPECT_EQ( square5[9]->qubit() , 2 ) ;  // 10

    EXPECT_NEAR( square5[0]->theta() , 1./ 1 , eps ) ;  // 01
    EXPECT_NEAR( square5[1]->theta() , 1./ 2 , eps ) ;  // 02

    EXPECT_NEAR( square5[2]->theta() , 1./ 3 , eps ) ;  // 03
    EXPECT_NEAR( square5[3]->theta() , 1./ 4 , eps ) ;  // 04

    EXPECT_NEAR( square5[4]->theta() , 1./ 5 , eps ) ;  // 05
    EXPECT_NEAR( square5[5]->theta() , 1./ 6 , eps ) ;  // 06

    EXPECT_NEAR( square5[6]->theta() , 1./ 7 , eps ) ;  // 07
    EXPECT_NEAR( square5[7]->theta() , 1./ 8 , eps ) ;  // 08

    EXPECT_NEAR( square5[8]->theta() , 1./ 9 , eps ) ;  // 09
    EXPECT_NEAR( square5[9]->theta() , 1./10 , eps ) ;  // 10

    auto square6 = test_f3c_SquareCircuit_init< T >( 6 ) ;

    EXPECT_EQ( square6[ 0]->qubit() , 0 ) ;  // 01
    EXPECT_EQ( square6[ 1]->qubit() , 2 ) ;  // 02
    EXPECT_EQ( square6[ 2]->qubit() , 4 ) ;  // 03

    EXPECT_EQ( square6[ 3]->qubit() , 1 ) ;  // 04
    EXPECT_EQ( square6[ 4]->qubit() , 3 ) ;  // 05

    EXPECT_EQ( square6[ 5]->qubit() , 0 ) ;  // 06
    EXPECT_EQ( square6[ 6]->qubit() , 2 ) ;  // 07
    EXPECT_EQ( square6[ 7]->qubit() , 4 ) ;  // 08

    EXPECT_EQ( square6[ 8]->qubit() , 1 ) ;  // 09
    EXPECT_EQ( square6[ 9]->qubit() , 3 ) ;  // 10

    EXPECT_EQ( square6[10]->qubit() , 0 ) ;  // 11
    EXPECT_EQ( square6[11]->qubit() , 2 ) ;  // 12
    EXPECT_EQ( square6[12]->qubit() , 4 ) ;  // 13

    EXPECT_EQ( square6[13]->qubit() , 1 ) ;  // 14
    EXPECT_EQ( square6[14]->qubit() , 3 ) ;  // 15

    EXPECT_NEAR( square6[ 0]->theta() , 1./ 1 , eps ) ;  // 01
    EXPECT_NEAR( square6[ 1]->theta() , 1./ 2 , eps ) ;  // 02
    EXPECT_NEAR( square6[ 2]->theta() , 1./ 3 , eps ) ;  // 03

    EXPECT_NEAR( square6[ 3]->theta() , 1./ 4 , eps ) ;  // 04
    EXPECT_NEAR( square6[ 4]->theta() , 1./ 5 , eps ) ;  // 05

    EXPECT_NEAR( square6[ 5]->theta() , 1./ 6 , eps ) ;  // 06
    EXPECT_NEAR( square6[ 6]->theta() , 1./ 7 , eps ) ;  // 07
    EXPECT_NEAR( square6[ 7]->theta() , 1./ 8 , eps ) ;  // 08

    EXPECT_NEAR( square6[ 8]->theta() , 1./ 9 , eps ) ;  // 09
    EXPECT_NEAR( square6[ 9]->theta() , 1./10 , eps ) ;  // 10

    EXPECT_NEAR( square6[10]->theta() , 1./11 , eps ) ;  // 11
    EXPECT_NEAR( square6[11]->theta() , 1./12 , eps ) ;  // 12
    EXPECT_NEAR( square6[12]->theta() , 1./13 , eps ) ;  // 13

    EXPECT_NEAR( square6[13]->theta() , 1./14 , eps ) ;  // 14
    EXPECT_NEAR( square6[14]->theta() , 1./15 , eps ) ;  // 15
  }


  //
  // firstIdx & lastIdx
  //
  {
    f3c::SquareCircuit< T >  square( 3 ) ;

    EXPECT_EQ( square.firstIdx( 0 ) , 0 ) ;
    EXPECT_EQ( square.firstIdx( 1 ) , 1 ) ;
    EXPECT_EQ( square.firstIdx( 2 ) , 2 ) ;

    EXPECT_EQ( square.lastIdx( 0 ) , 0 ) ;
    EXPECT_EQ( square.lastIdx( 1 ) , 1 ) ;
    EXPECT_EQ( square.lastIdx( 2 ) , 2 ) ;
  }

  {
    f3c::SquareCircuit< T >  square( 4 ) ;

    EXPECT_EQ( square.firstIdx( 0 ) , 0 ) ;
    EXPECT_EQ( square.firstIdx( 1 ) , 2 ) ;
    EXPECT_EQ( square.firstIdx( 2 ) , 3 ) ;
    EXPECT_EQ( square.firstIdx( 3 ) , 5 ) ;

    EXPECT_EQ( square.lastIdx( 0 ) , 1 ) ;
    EXPECT_EQ( square.lastIdx( 1 ) , 2 ) ;
    EXPECT_EQ( square.lastIdx( 2 ) , 4 ) ;
    EXPECT_EQ( square.lastIdx( 3 ) , 5 ) ;
  }

  {
    f3c::SquareCircuit< T >  square( 5 ) ;

    EXPECT_EQ( square.firstIdx( 0 ) , 0 ) ;
    EXPECT_EQ( square.firstIdx( 1 ) , 2 ) ;
    EXPECT_EQ( square.firstIdx( 2 ) , 4 ) ;
    EXPECT_EQ( square.firstIdx( 3 ) , 6 ) ;
    EXPECT_EQ( square.firstIdx( 4 ) , 8 ) ;

    EXPECT_EQ( square.lastIdx( 0 ) , 1 ) ;
    EXPECT_EQ( square.lastIdx( 1 ) , 3 ) ;
    EXPECT_EQ( square.lastIdx( 2 ) , 5 ) ;
    EXPECT_EQ( square.lastIdx( 3 ) , 7 ) ;
    EXPECT_EQ( square.lastIdx( 4 ) , 9 ) ;
  }

  {
    f3c::SquareCircuit< T >  square( 6 ) ;

    EXPECT_EQ( square.firstIdx( 0 ) ,  0 ) ;
    EXPECT_EQ( square.firstIdx( 1 ) ,  3 ) ;
    EXPECT_EQ( square.firstIdx( 2 ) ,  5 ) ;
    EXPECT_EQ( square.firstIdx( 3 ) ,  8 ) ;
    EXPECT_EQ( square.firstIdx( 4 ) , 10 ) ;
    EXPECT_EQ( square.firstIdx( 5 ) , 13 ) ;

    EXPECT_EQ( square.lastIdx( 0 ) ,  2 ) ;
    EXPECT_EQ( square.lastIdx( 1 ) ,  4 ) ;
    EXPECT_EQ( square.lastIdx( 2 ) ,  7 ) ;
    EXPECT_EQ( square.lastIdx( 3 ) ,  9 ) ;
    EXPECT_EQ( square.lastIdx( 4 ) , 12 ) ;
    EXPECT_EQ( square.lastIdx( 5 ) , 14 ) ;
  }

  {
    f3c::SquareCircuit< T >  square( 7 ) ;

    EXPECT_EQ( square.firstIdx( 0 ) ,  0 ) ;
    EXPECT_EQ( square.firstIdx( 1 ) ,  3 ) ;
    EXPECT_EQ( square.firstIdx( 2 ) ,  6 ) ;
    EXPECT_EQ( square.firstIdx( 3 ) ,  9 ) ;
    EXPECT_EQ( square.firstIdx( 4 ) , 12 ) ;
    EXPECT_EQ( square.firstIdx( 5 ) , 15 ) ;
    EXPECT_EQ( square.firstIdx( 6 ) , 18 ) ;

    EXPECT_EQ( square.lastIdx( 0 ) ,  2 ) ;
    EXPECT_EQ( square.lastIdx( 1 ) ,  5 ) ;
    EXPECT_EQ( square.lastIdx( 2 ) ,  8 ) ;
    EXPECT_EQ( square.lastIdx( 3 ) , 11 ) ;
    EXPECT_EQ( square.lastIdx( 4 ) , 14 ) ;
    EXPECT_EQ( square.lastIdx( 5 ) , 17 ) ;
    EXPECT_EQ( square.lastIdx( 6 ) , 20 ) ;
  }

  {
    f3c::SquareCircuit< T >  square( 8 ) ;

    EXPECT_EQ( square.firstIdx( 0 ) ,  0 ) ;
    EXPECT_EQ( square.firstIdx( 1 ) ,  4 ) ;
    EXPECT_EQ( square.firstIdx( 2 ) ,  7 ) ;
    EXPECT_EQ( square.firstIdx( 3 ) , 11 ) ;
    EXPECT_EQ( square.firstIdx( 4 ) , 14 ) ;
    EXPECT_EQ( square.firstIdx( 5 ) , 18 ) ;
    EXPECT_EQ( square.firstIdx( 6 ) , 21 ) ;
    EXPECT_EQ( square.firstIdx( 7 ) , 25 ) ;

    EXPECT_EQ( square.lastIdx( 0 ) ,  3 ) ;
    EXPECT_EQ( square.lastIdx( 1 ) ,  6 ) ;
    EXPECT_EQ( square.lastIdx( 2 ) , 10 ) ;
    EXPECT_EQ( square.lastIdx( 3 ) , 13 ) ;
    EXPECT_EQ( square.lastIdx( 4 ) , 17 ) ;
    EXPECT_EQ( square.lastIdx( 5 ) , 20 ) ;
    EXPECT_EQ( square.lastIdx( 6 ) , 24 ) ;
    EXPECT_EQ( square.lastIdx( 7 ) , 27 ) ;
  }


  //
  // toTriangle
  //
  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 3 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 3 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 8*eps ) ;
  }

  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 4 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 4 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 16*eps ) ;
  }

  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 5 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 5 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 32*eps ) ;
  }

  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 6 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 6 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 64*eps ) ;
  }

  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 7 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 7 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 128*eps ) ;
  }

  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 8 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 8 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 256*eps ) ;
  }

  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 9 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 9 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 512*eps ) ;
  }

  {
    auto square = test_f3c_SquareCircuit_initXY< T >( 10 ) ;
    auto check  = test_f3c_SquareCircuit_initXY< T >( 10 ) ;

    auto triangle = square.toTriangle() ;
    EXPECT_NEAR( qclab::nrmF( triangle , check ) , 0.0 , 1024*eps ) ;
  }

}


/*
 * float
 */
TEST( f3c_SquareCircuit , complex_float ) {
  test_f3c_SquareCircuit< std::complex< float > >() ;
}

/*
 * double
 */
TEST( f3c_SquareCircuit , complex_double ) {
  test_f3c_SquareCircuit< std::complex< double > >() ;
}

