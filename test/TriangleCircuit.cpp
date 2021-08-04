#include <gtest/gtest.h>
#include "f3c/TriangleCircuit.hpp"
#include "f3c/SquareCircuit.hpp"
#include "qclab/qgates/RotationXX.hpp"
#include "f3c/qgates/RotationXY.hpp"

template <typename T>
auto test_f3c_TriangleCircuit_init( const int n ) {

  using XX = qclab::qgates::RotationXX< T > ;

  f3c::TriangleCircuit< T , XX >  triangle( n ) ;
  int c = 0 ;
  for ( int l = 0; l < n-1; l++ ) {
    for ( int i = 0; i < n-l-1; i++ ) {
      triangle[c] = std::make_unique< XX >( n-i-2 , n-i-1 , 1./(c+1) ) ;
      c++ ;
    }
  }
  return triangle ;

}


template <typename T>
auto test_f3c_TriangleCircuit_initXY( const int n ) {

  using XY = f3c::qgates::RotationXY< T > ;

  f3c::TriangleCircuit< T , XY >  triangle( n ) ;
  int c = 0 ;
  for ( int l = 0; l < n-1; l++ ) {
    for ( int i = 0; i < n-l-1; i++ ) {
      triangle[c] = std::make_unique< XY >( n-i-2 , n-i-1 ,
                                            1./(c+1) , 1. + 1./(c+1) ) ;
      c++ ;
    }
  }
  return triangle ;

}


template <typename T>
void test_f3c_TriangleCircuit() {

  using R = qclab::real_t< T > ;
  using G2 = qclab::qgates::QGate2< T > ;
  using XX = qclab::qgates::RotationXX< T > ;
  using XY =   f3c::qgates::RotationXY< T > ;

  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    f3c::TriangleCircuit< T >  triangle( 3 ) ;
    EXPECT_EQ( triangle.nbQubits() , 3 ) ;
    EXPECT_EQ( triangle.capacity() , 3 ) ;
    EXPECT_TRUE( triangle.ascend() ) ;
    EXPECT_FALSE( triangle.descend() ) ;
  }

  {
    f3c::TriangleCircuit< T >  triangle2( 2 ) ;
    EXPECT_EQ( triangle2.nbQubits() , 2 ) ;
    EXPECT_EQ( triangle2.capacity() , 1 ) ;

    f3c::TriangleCircuit< T >  triangle3( 3 ) ;
    EXPECT_EQ( triangle3.nbQubits() , 3 ) ;
    EXPECT_EQ( triangle3.capacity() , 3 ) ;

    f3c::TriangleCircuit< T >  triangle4( 4 ) ;
    EXPECT_EQ( triangle4.nbQubits() , 4 ) ;
    EXPECT_EQ( triangle4.capacity() , 6 ) ;

    f3c::TriangleCircuit< T >  triangle5( 5 ) ;
    EXPECT_EQ( triangle5.nbQubits() , 5 ) ;
    EXPECT_EQ( triangle5.capacity() , 10 ) ;

    f3c::TriangleCircuit< T >  triangle6( 6 ) ;
    EXPECT_EQ( triangle6.nbQubits() , 6 ) ;
    EXPECT_EQ( triangle6.capacity() , 15 ) ;

    f3c::TriangleCircuit< T >  triangle7( 7 ) ;
    EXPECT_EQ( triangle7.nbQubits() , 7 ) ;
    EXPECT_EQ( triangle7.capacity() , 21 ) ;

    f3c::TriangleCircuit< T >  triangle8( 8 ) ;
    EXPECT_EQ( triangle8.nbQubits() , 8 ) ;
    EXPECT_EQ( triangle8.capacity() , 28 ) ;
  }

  {
    int n = 6 ;
    f3c::TriangleCircuit< T , XX >  triangle( n ) ;
    int c = 0 ;
    for ( int l = 0; l < n-1; l++ ) {
      for ( int i = 0; i < n-l-1; i++ ) {
        triangle[c] = std::make_unique< XX >( n-i-2 , n-i-1 , 1./(c+1) ) ;
        c++ ;
      }
    }
    auto matA = triangle.matrix() ;

    // ascending
    EXPECT_TRUE(  triangle.ascend() ) ;
    EXPECT_FALSE( triangle.descend() ) ;

    EXPECT_EQ( triangle[ 0]->qubit() , 4 ) ;  // 01
    EXPECT_EQ( triangle[ 1]->qubit() , 3 ) ;  // 02
    EXPECT_EQ( triangle[ 2]->qubit() , 2 ) ;  // 03
    EXPECT_EQ( triangle[ 3]->qubit() , 1 ) ;  // 04
    EXPECT_EQ( triangle[ 4]->qubit() , 0 ) ;  // 05

    EXPECT_EQ( triangle[ 5]->qubit() , 4 ) ;  // 06
    EXPECT_EQ( triangle[ 6]->qubit() , 3 ) ;  // 07
    EXPECT_EQ( triangle[ 7]->qubit() , 2 ) ;  // 08
    EXPECT_EQ( triangle[ 8]->qubit() , 1 ) ;  // 09

    EXPECT_EQ( triangle[ 9]->qubit() , 4 ) ;  // 10
    EXPECT_EQ( triangle[10]->qubit() , 3 ) ;  // 11
    EXPECT_EQ( triangle[11]->qubit() , 2 ) ;  // 12

    EXPECT_EQ( triangle[12]->qubit() , 4 ) ;  // 13
    EXPECT_EQ( triangle[13]->qubit() , 3 ) ;  // 14

    EXPECT_EQ( triangle[14]->qubit() , 4 ) ;  // 15

    for ( int i = 0; i < triangle.nbGates(); i++ ) {
      const int c = i + 1 ;
      EXPECT_NEAR( triangle[i]->theta() , 1./c , 10*eps ) ;
    }

    // makeDescend
    triangle.makeDescend() ;
    EXPECT_FALSE( triangle.ascend() ) ;
    EXPECT_TRUE(  triangle.descend() ) ;

    EXPECT_EQ( triangle[ 0]->qubit() , 4 ) ;  // 01

    EXPECT_EQ( triangle[ 1]->qubit() , 3 ) ;  // 02
    EXPECT_EQ( triangle[ 2]->qubit() , 4 ) ;  // 03

    EXPECT_EQ( triangle[ 3]->qubit() , 2 ) ;  // 04
    EXPECT_EQ( triangle[ 4]->qubit() , 3 ) ;  // 05
    EXPECT_EQ( triangle[ 5]->qubit() , 4 ) ;  // 06

    EXPECT_EQ( triangle[ 6]->qubit() , 1 ) ;  // 07
    EXPECT_EQ( triangle[ 7]->qubit() , 2 ) ;  // 08
    EXPECT_EQ( triangle[ 8]->qubit() , 3 ) ;  // 09
    EXPECT_EQ( triangle[ 9]->qubit() , 4 ) ;  // 10

    EXPECT_EQ( triangle[10]->qubit() , 0 ) ;  // 11
    EXPECT_EQ( triangle[11]->qubit() , 1 ) ;  // 12
    EXPECT_EQ( triangle[12]->qubit() , 2 ) ;  // 13
    EXPECT_EQ( triangle[13]->qubit() , 3 ) ;  // 14
    EXPECT_EQ( triangle[14]->qubit() , 4 ) ;  // 15

    EXPECT_NEAR( triangle[ 0]->theta() , 1./ 1 , 10*eps ) ;  // 01

    EXPECT_NEAR( triangle[ 1]->theta() , 1./ 2 , 10*eps ) ;  // 02
    EXPECT_NEAR( triangle[ 2]->theta() , 1./ 6 , 10*eps ) ;  // 03

    EXPECT_NEAR( triangle[ 3]->theta() , 1./ 3 , 10*eps ) ;  // 04
    EXPECT_NEAR( triangle[ 4]->theta() , 1./ 7 , 10*eps ) ;  // 05
    EXPECT_NEAR( triangle[ 5]->theta() , 1./10 , 10*eps ) ;  // 06

    EXPECT_NEAR( triangle[ 6]->theta() , 1./ 4 , 10*eps ) ;  // 07
    EXPECT_NEAR( triangle[ 7]->theta() , 1./ 8 , 10*eps ) ;  // 08
    EXPECT_NEAR( triangle[ 8]->theta() , 1./11 , 10*eps ) ;  // 09
    EXPECT_NEAR( triangle[ 9]->theta() , 1./13 , 10*eps ) ;  // 10

    EXPECT_NEAR( triangle[10]->theta() , 1./ 5 , 10*eps ) ;  // 11
    EXPECT_NEAR( triangle[11]->theta() , 1./ 9 , 10*eps ) ;  // 12
    EXPECT_NEAR( triangle[12]->theta() , 1./12 , 10*eps ) ;  // 13
    EXPECT_NEAR( triangle[13]->theta() , 1./14 , 10*eps ) ;  // 14
    EXPECT_NEAR( triangle[14]->theta() , 1./15 , 10*eps ) ;  // 15

    auto matD = triangle.matrix() ;
    for ( int c = 0; c < matA.cols(); c++ ) {
      for ( int r = 0; r < matA.rows(); r++ ) {
        EXPECT_NEAR( std::real( matA(r,c) ) , std::real( matD(r,c) ) , 10*eps );
        EXPECT_NEAR( std::imag( matA(r,c) ) , std::imag( matD(r,c) ) , 10*eps );
      }
    }

    // makeAscend
    triangle.makeAscend() ;
    EXPECT_TRUE(  triangle.ascend() ) ;
    EXPECT_FALSE( triangle.descend() ) ;

    EXPECT_EQ( triangle[ 0]->qubit() , 4 ) ;  // 01
    EXPECT_EQ( triangle[ 1]->qubit() , 3 ) ;  // 02
    EXPECT_EQ( triangle[ 2]->qubit() , 2 ) ;  // 03
    EXPECT_EQ( triangle[ 3]->qubit() , 1 ) ;  // 04
    EXPECT_EQ( triangle[ 4]->qubit() , 0 ) ;  // 05

    EXPECT_EQ( triangle[ 5]->qubit() , 4 ) ;  // 06
    EXPECT_EQ( triangle[ 6]->qubit() , 3 ) ;  // 07
    EXPECT_EQ( triangle[ 7]->qubit() , 2 ) ;  // 08
    EXPECT_EQ( triangle[ 8]->qubit() , 1 ) ;  // 09

    EXPECT_EQ( triangle[ 9]->qubit() , 4 ) ;  // 10
    EXPECT_EQ( triangle[10]->qubit() , 3 ) ;  // 11
    EXPECT_EQ( triangle[11]->qubit() , 2 ) ;  // 12

    EXPECT_EQ( triangle[12]->qubit() , 4 ) ;  // 13
    EXPECT_EQ( triangle[13]->qubit() , 3 ) ;  // 14

    EXPECT_EQ( triangle[14]->qubit() , 4 ) ;  // 15

    EXPECT_NEAR( triangle[ 0]->theta() , 1./ 1 , 10*eps ) ;  // 01

    EXPECT_NEAR( triangle[ 1]->theta() , 1./ 2 , 10*eps ) ;  // 02
    EXPECT_NEAR( triangle[ 2]->theta() , 1./ 3 , 10*eps ) ;  // 03

    EXPECT_NEAR( triangle[ 3]->theta() , 1./ 4 , 10*eps ) ;  // 04
    EXPECT_NEAR( triangle[ 4]->theta() , 1./ 5 , 10*eps ) ;  // 05
    EXPECT_NEAR( triangle[ 5]->theta() , 1./ 6 , 10*eps ) ;  // 06

    EXPECT_NEAR( triangle[ 6]->theta() , 1./ 7 , 10*eps ) ;  // 07
    EXPECT_NEAR( triangle[ 7]->theta() , 1./ 8 , 10*eps ) ;  // 08
    EXPECT_NEAR( triangle[ 8]->theta() , 1./ 9 , 10*eps ) ;  // 09
    EXPECT_NEAR( triangle[ 9]->theta() , 1./10 , 10*eps ) ;  // 10

    EXPECT_NEAR( triangle[10]->theta() , 1./11 , 10*eps ) ;  // 11
    EXPECT_NEAR( triangle[11]->theta() , 1./12 , 10*eps ) ;  // 12
    EXPECT_NEAR( triangle[12]->theta() , 1./13 , 10*eps ) ;  // 13
    EXPECT_NEAR( triangle[13]->theta() , 1./14 , 10*eps ) ;  // 14
    EXPECT_NEAR( triangle[14]->theta() , 1./15 , 10*eps ) ;  // 15

    auto matA2 = triangle.matrix() ;
    for ( int c = 0; c < matA.cols(); c++ ) {
      for ( int r = 0; r < matA.rows(); r++ ) {
        EXPECT_EQ( std::real( matA(r,c) ) , std::real( matA2(r,c) ) ) ;
        EXPECT_EQ( std::imag( matA(r,c) ) , std::imag( matA2(r,c) ) ) ;
      }
    }
  }

  {
    auto tr2a = test_f3c_TriangleCircuit_init< T >( 2 ) ;
    auto tr2b = test_f3c_TriangleCircuit_init< T >( 2 ) ; tr2b.makeDescend() ;
    EXPECT_NEAR( qclab::nrmF( tr2a , tr2b ) , 0.0 , 10*eps ) ;
    tr2b.makeAscend() ;
    EXPECT_NEAR( qclab::nrmF( tr2a , tr2b ) , 0.0 , 0 ) ;

    auto tr3a = test_f3c_TriangleCircuit_init< T >( 3 ) ;
    auto tr3b = test_f3c_TriangleCircuit_init< T >( 3 ) ; tr3b.makeDescend() ;
    EXPECT_NEAR( qclab::nrmF( tr3a , tr3b ) , 0.0 , 10*eps ) ;
    tr3b.makeAscend() ;
    EXPECT_NEAR( qclab::nrmF( tr3a , tr3b ) , 0.0 , 0 ) ;

    auto tr7a = test_f3c_TriangleCircuit_init< T >( 7 ) ;
    auto tr7b = test_f3c_TriangleCircuit_init< T >( 7 ) ; tr7b.makeDescend() ;
    EXPECT_NEAR( qclab::nrmF( tr7a , tr7b ) , 0.0 , 100*eps ) ;
    tr7b.makeAscend() ;
    EXPECT_NEAR( qclab::nrmF( tr7a , tr7b ) , 0.0 , 0 ) ;

    auto tr8a = test_f3c_TriangleCircuit_init< T >( 8 ) ;
    auto tr8b = test_f3c_TriangleCircuit_init< T >( 8 ) ; tr8b.makeDescend() ;
    EXPECT_NEAR( qclab::nrmF( tr8a , tr8b ) , 0.0 , 100*eps ) ;
    tr8b.makeAscend() ;
    EXPECT_NEAR( qclab::nrmF( tr8a , tr8b ) , 0.0 , 0 ) ;
  }

  {
    int n = 6 ;
    auto triangle = test_f3c_TriangleCircuit_init< T >( n ) ;

    // ascIdx
    EXPECT_EQ( triangle.ascIdx( 0 , 0 ) ,  4 ) ;
    EXPECT_EQ( triangle.ascIdx( 0 , 1 ) ,  3 ) ;
    EXPECT_EQ( triangle.ascIdx( 0 , 2 ) ,  2 ) ;
    EXPECT_EQ( triangle.ascIdx( 0 , 3 ) ,  1 ) ;
    EXPECT_EQ( triangle.ascIdx( 0 , 4 ) ,  0 ) ;

    EXPECT_EQ( triangle.ascIdx( 1 , 1 ) ,  8 ) ;
    EXPECT_EQ( triangle.ascIdx( 1 , 2 ) ,  7 ) ;
    EXPECT_EQ( triangle.ascIdx( 1 , 3 ) ,  6 ) ;
    EXPECT_EQ( triangle.ascIdx( 1 , 4 ) ,  5 ) ;

    EXPECT_EQ( triangle.ascIdx( 2 , 2 ) , 11 ) ;
    EXPECT_EQ( triangle.ascIdx( 2 , 3 ) , 10 ) ;
    EXPECT_EQ( triangle.ascIdx( 2 , 4 ) ,  9 ) ;

    EXPECT_EQ( triangle.ascIdx( 3 , 3 ) , 13 ) ;
    EXPECT_EQ( triangle.ascIdx( 3 , 4 ) , 12 ) ;

    EXPECT_EQ( triangle.ascIdx( 4 , 4 ) , 14 ) ;

    //std::cout << std::endl ;
    //for ( int l = 0; l < n-1; l++ ) {
    //  for ( int q = l; q < n-1; q++ ) {
    //    std::cout << "* l = " << l << ", q = " << q << ": "
    //              << triangle.ascIdx( l , q ) << std::endl ;
    //  }
    //  std::cout << std::endl ;
    //}

    triangle.makeDescend() ;

    // desIdx
    EXPECT_EQ( triangle.desIdx( 0 , 4 ) ,  0 ) ;

    EXPECT_EQ( triangle.desIdx( 1 , 3 ) ,  1 ) ;
    EXPECT_EQ( triangle.desIdx( 1 , 4 ) ,  2 ) ;

    EXPECT_EQ( triangle.desIdx( 2 , 2 ) ,  3 ) ;
    EXPECT_EQ( triangle.desIdx( 2 , 3 ) ,  4 ) ;
    EXPECT_EQ( triangle.desIdx( 2 , 4 ) ,  5 ) ;

    EXPECT_EQ( triangle.desIdx( 3 , 1 ) ,  6 ) ;
    EXPECT_EQ( triangle.desIdx( 3 , 2 ) ,  7 ) ;
    EXPECT_EQ( triangle.desIdx( 3 , 3 ) ,  8 ) ;
    EXPECT_EQ( triangle.desIdx( 3 , 4 ) ,  9 ) ;

    EXPECT_EQ( triangle.desIdx( 4 , 0 ) , 10 ) ;
    EXPECT_EQ( triangle.desIdx( 4 , 1 ) , 11 ) ;
    EXPECT_EQ( triangle.desIdx( 4 , 2 ) , 12 ) ;
    EXPECT_EQ( triangle.desIdx( 4 , 3 ) , 13 ) ;
    EXPECT_EQ( triangle.desIdx( 4 , 4 ) , 14 ) ;

    //std::cout << std::endl ;
    //for ( int l = 0; l < n-1; l++ ) {
    //  for ( int q = n-l-2; q < n-1; q++ ) {
    //    std::cout << "* l = " << l << ", q = " << q << ": "
    //              << triangle.desIdx( l , q ) << std::endl ;
    //  }
    //  std::cout << std::endl ;
    //}
  }


  //
  // Merge (ascending + left)
  //
  {
    int n = 6 ;
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( n ) ;
    qclab::QCircuit< T , qclab::qgates::QGate2< T > >  circheck( n ) ;
    for ( auto it = triangle.begin(); it != triangle.end(); ++it ) {
      circheck.push_back( std::make_unique< XY >( **it ) ) ;
    }
    EXPECT_EQ( qclab::nrmF( triangle , circheck ) , 0.0 ) ;

    // merge gate at qubits 0 and 1
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 0 , 1 , 1.0 , 1.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    XY gate( 0 , 1 , 1.0 , 1.5 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 1 and 2
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 1 , 2 , 0.5 , 1.1 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    int qnew[2] = { 1 , 2 } ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.5 , 1.1 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 2 and 3
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 2 , 3 , -.4 , -.3 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 2 ; qnew[1] = 3 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( -.4 , -.3 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 3 and 4
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 3 , 4 , 0.1 , 0.2 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 3 ; qnew[1] = 4 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.1 , 0.2 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 4 and 5
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 4 , 5 , 0.3 , 0.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 4 ; qnew[1] = 5 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.3 , 0.5 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;
  }


  //
  // Merge (ascending + right)
  //
  {
    int n = 6 ;
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( n ) ;
    qclab::QCircuit< T , qclab::qgates::QGate2< T > >  circheck( n ) ;
    for ( auto it = triangle.begin(); it != triangle.end(); ++it ) {
      circheck.push_back( std::make_unique< XY >( **it ) ) ;
    }
    EXPECT_EQ( qclab::nrmF( triangle , circheck ) , 0.0 ) ;

    // merge gate at qubits 0 and 1
    circheck.push_back( std::make_unique< XY >( 0 , 1 , 1.0 , 1.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    XY gate( 0 , 1 , 1.0 , 1.5 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 1 and 2
    circheck.push_back( std::make_unique< XY >( 1 , 2 , 0.5 , 1.1 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    int qnew[2] = { 1 , 2 } ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.5 , 1.1 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 2 and 3
    circheck.push_back( std::make_unique< XY >( 2 , 3 , -.4 , -.3 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 2 ; qnew[1] = 3 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( -.4 , -.3 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 3 and 4
    circheck.push_back( std::make_unique< XY >( 3 , 4 , 0.1 , 0.2 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 3 ; qnew[1] = 4 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.1 , 0.2 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 4 and 5
    circheck.push_back( std::make_unique< XY >( 4 , 5 , 0.3 , 0.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 4 ; qnew[1] = 5 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.3 , 0.5 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;
  }


  //
  // Merge (descending + left)
  //
  {
    int n = 6 ;
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( n ) ;
    triangle.makeDescend() ;
    qclab::QCircuit< T , qclab::qgates::QGate2< T > >  circheck( n ) ;
    for ( auto it = triangle.begin(); it != triangle.end(); ++it ) {
      circheck.push_back( std::make_unique< XY >( **it ) ) ;
    }
    EXPECT_EQ( qclab::nrmF( triangle , circheck ) , 0.0 ) ;

    // merge gate at qubits 0 and 1
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 0 , 1 , 1.0 , 1.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    XY gate( 0 , 1 , 1.0 , 1.5 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 1 and 2
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 1 , 2 , 0.5 , 1.1 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    int qnew[2] = { 1 , 2 } ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.5 , 1.1 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 2 and 3
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 2 , 3 , -.4 , -.3 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 2 ; qnew[1] = 3 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( -.4 , -.3 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 3 and 4
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 3 , 4 , 0.1 , 0.2 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 3 ; qnew[1] = 4 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.1 , 0.2 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 4 and 5
    circheck.insert( circheck.begin() ,
                     std::make_unique< XY >( 4 , 5 , 0.3 , 0.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 4 ; qnew[1] = 5 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.3 , 0.5 ) ;
    triangle.merge( qclab::Side::Left , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;
  }


  //
  // Merge (descending + right)
  //
  {
    int n = 6 ;
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( n ) ;
    triangle.makeDescend() ;
    qclab::QCircuit< T , qclab::qgates::QGate2< T > >  circheck( n ) ;
    for ( auto it = triangle.begin(); it != triangle.end(); ++it ) {
      circheck.push_back( std::make_unique< XY >( **it ) ) ;
    }
    EXPECT_EQ( qclab::nrmF( triangle , circheck ) , 0.0 ) ;

    // merge gate at qubits 0 and 1
    circheck.push_back( std::make_unique< XY >( 0 , 1 , 1.0 , 1.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    XY gate( 0 , 1 , 1.0 , 1.5 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 1 and 2
    circheck.push_back( std::make_unique< XY >( 1 , 2 , 0.5 , 1.1 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    int qnew[2] = { 1 , 2 } ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.5 , 1.1 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 2 and 3
    circheck.push_back( std::make_unique< XY >( 2 , 3 , -.4 , -.3 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 2 ; qnew[1] = 3 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( -.4 , -.3 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 3 and 4
    circheck.push_back( std::make_unique< XY >( 3 , 4 , 0.1 , 0.2 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 3 ; qnew[1] = 4 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.1 , 0.2 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;

    // merge gate at qubits 4 and 5
    circheck.push_back( std::make_unique< XY >( 4 , 5 , 0.3 , 0.5 ) ) ;
    EXPECT_TRUE( qclab::nrmF( triangle , circheck ) > 100*eps ) ;
    qnew[0] = 4 ; qnew[1] = 5 ;
    gate.setQubits( &qnew[0] ) ;
    gate.update( 0.3 , 0.5 ) ;
    triangle.merge( qclab::Side::Right , gate ) ;
    EXPECT_NEAR( qclab::nrmF( triangle , circheck ) , 0.0 , 100*eps ) ;
  }


  //
  // toSquare
  //
  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 3 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 3 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 8*10*eps ) ;
  }

  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 4 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 4 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 16*10*eps ) ;
  }

  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 5 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 5 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 32*10*eps ) ;
  }

  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 6 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 6 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 64*10*eps ) ;
  }

  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 7 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 7 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 128*10*eps ) ;
  }

  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 8 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 8 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 256*10*eps ) ;
  }

  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 9 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 9 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 512*10*eps ) ;
  }

  {
    auto triangle = test_f3c_TriangleCircuit_initXY< T >( 10 ) ;
    auto check    = test_f3c_TriangleCircuit_initXY< T >( 10 ) ;

    auto square = triangle.toSquare() ;
    EXPECT_NEAR( qclab::nrmF( square , check ) , 0.0 , 1024*10*eps ) ;
  }

}


/*
 * float
 */
TEST( f3c_TriangleCircuit , complex_float ) {
  test_f3c_TriangleCircuit< std::complex< float > >() ;
}

/*
 * double
 */
TEST( f3c_TriangleCircuit , complex_double ) {
  test_f3c_TriangleCircuit< std::complex< double > >() ;
}

