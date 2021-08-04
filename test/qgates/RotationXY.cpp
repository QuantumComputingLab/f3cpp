#include <gtest/gtest.h>
#include "f3c/qgates/RotationXY.hpp"

template <typename T>
void test_f3c_qgates_RotationXY() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    f3c::qgates::RotationXY< T >  Rxy ;

    EXPECT_EQ( Rxy.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_FALSE( Rxy.fixed() ) ;       // fixed
    EXPECT_FALSE( Rxy.controlled() ) ;  // controlled

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( Rxy.matrix() == eye ) ;

    // qubits
    EXPECT_EQ( Rxy.qubits().size() , 2 ) ;
    EXPECT_EQ( Rxy.qubits()[0] , 0 ) ;
    EXPECT_EQ( Rxy.qubits()[1] , 1 ) ;
    int qnew[2] = { 3 , 5 } ;
    Rxy.setQubits( &qnew[0] ) ;
    EXPECT_EQ( Rxy.qubits().size() , 2 ) ;
    EXPECT_EQ( Rxy.qubits()[0] , 3 ) ;
    EXPECT_EQ( Rxy.qubits()[1] , 5 ) ;

    // rotations and thetas
    auto [ rot0 , rot1 ] = Rxy.rotations() ;
    EXPECT_EQ( rot0.theta() , 0 ) ;
    EXPECT_EQ( rot1.theta() , 0 ) ;
    auto [ theta0 , theta1 ] = Rxy.thetas() ;
    EXPECT_EQ( theta0 , 0 ) ;
    EXPECT_EQ( theta1 , 0 ) ;

    // update(rotation0,rotation1)
    qclab::QRotation< R >  new_rot0( 0.5 ) ;
    qclab::QRotation< R >  new_rot1( 0.3 ) ;
    Rxy.update( new_rot0 , new_rot1 ) ;
    {
      auto [ r0 , r1 ] = Rxy.rotations() ;
      EXPECT_TRUE( r0 == new_rot0 ) ;
      EXPECT_TRUE( r1 == new_rot1 ) ;
      auto [ theta0 , theta1 ] = Rxy.thetas() ;
      EXPECT_NEAR( theta0 , 0.5 , eps ) ;
      EXPECT_NEAR( theta1 , 0.3 , eps ) ;
    }

    // update(theta0,theta1)
    theta0 = pi/2 ;
    theta1 = pi/3 ;
    Rxy.update( theta0 , theta1 ) ;
    {
      auto [ theta0 , theta1 ] = Rxy.thetas() ;
      EXPECT_NEAR( theta0 , theta0 , eps ) ;
      EXPECT_NEAR( theta1 , theta1 , eps ) ;
      auto [ r0 , r1 ] = Rxy.rotations() ;
      EXPECT_NEAR( r0.theta() , theta0 , eps ) ;
      EXPECT_NEAR( r0.cos() , std::cos( theta0/2 ) , eps ) ;
      EXPECT_NEAR( r0.sin() , std::sin( theta0/2 ) , eps ) ;
      EXPECT_NEAR( r1.theta() , theta1 , eps ) ;
      EXPECT_NEAR( r1.cos() , std::cos( theta1/2 ) , eps ) ;
      EXPECT_NEAR( r1.sin() , std::sin( theta1/2 ) , eps ) ;
    }

    // matrix
    const R cos_plus  = std::cos( theta0/2 + theta1/2 ) ;
    const R sin_plus  = std::sin( theta0/2 + theta1/2 ) ;
    const R cos_minus = std::cos( theta0/2 - theta1/2 ) ;
    const R sin_minus = std::sin( theta0/2 - theta1/2 ) ;

    EXPECT_NEAR( std::real( Rxy.matrix()(0,0) ) , cos_minus , eps ) ;
    EXPECT_NEAR( std::real( Rxy.matrix()(1,1) ) , cos_plus  , eps ) ;
    EXPECT_NEAR( std::real( Rxy.matrix()(2,2) ) , cos_plus  , eps ) ;
    EXPECT_NEAR( std::real( Rxy.matrix()(3,3) ) , cos_minus , eps ) ;
    EXPECT_EQ( std::imag( Rxy.matrix()(0,0) ) , T(0) ) ;
    EXPECT_EQ( std::imag( Rxy.matrix()(1,1) ) , T(0) ) ;
    EXPECT_EQ( std::imag( Rxy.matrix()(2,2) ) , T(0) ) ;
    EXPECT_EQ( std::imag( Rxy.matrix()(3,3) ) , T(0) ) ;

    EXPECT_EQ( std::real( Rxy.matrix()(3,0) ) , T(0) ) ;
    EXPECT_EQ( std::real( Rxy.matrix()(2,1) ) , T(0) ) ;
    EXPECT_EQ( std::real( Rxy.matrix()(1,2) ) , T(0) ) ;
    EXPECT_EQ( std::real( Rxy.matrix()(0,3) ) , T(0) ) ;
    EXPECT_NEAR( std::imag( Rxy.matrix()(3,0) ) , -sin_minus , eps ) ;
    EXPECT_NEAR( std::imag( Rxy.matrix()(2,1) ) , -sin_plus  , eps ) ;
    EXPECT_NEAR( std::imag( Rxy.matrix()(1,2) ) , -sin_plus  , eps ) ;
    EXPECT_NEAR( std::imag( Rxy.matrix()(0,3) ) , -sin_minus , eps ) ;

    EXPECT_EQ( Rxy.matrix()(1,0) , T(0) ) ;
    EXPECT_EQ( Rxy.matrix()(2,0) , T(0) ) ;
    EXPECT_EQ( Rxy.matrix()(0,1) , T(0) ) ;
    EXPECT_EQ( Rxy.matrix()(3,1) , T(0) ) ;
    EXPECT_EQ( Rxy.matrix()(0,2) , T(0) ) ;
    EXPECT_EQ( Rxy.matrix()(3,2) , T(0) ) ;
    EXPECT_EQ( Rxy.matrix()(1,3) , T(0) ) ;
    EXPECT_EQ( Rxy.matrix()(2,3) , T(0) ) ;

    // print
    Rxy.print() ;

    // toQASM
    {
      std::stringstream qasm ;
      EXPECT_EQ( Rxy.toQASM( qasm ) , 0 ) ;
      auto str_pi = qclab::qasm( pi/2 ) ;
      auto [ theta0 , theta1 ] = Rxy.thetas() ;
      auto str_theta0 = qclab::qasm( theta0 ) ;
      auto str_theta1 = qclab::qasm( theta1 ) ;
      std::string qasm_check =
        "rx(" + str_pi + ") q[3];\n"     + "rx(" + str_pi + ") q[5];\n" +
        "cx q[3], q[5];\n" +
        "rx(" + str_theta0 + ") q[3];\n" + "rz(" + str_theta1 + ") q[5];\n" +
        "cx q[3], q[5];\n" +
        "rx(-" + str_pi + ") q[3];\n"    + "rx(-" + str_pi + ") q[5];\n" ;
      EXPECT_EQ( qasm.str() , qasm_check ) ;
      std::cout << qasm.str() ;
    }

    // operators == and !=
    {
      qclab::QRotation< R >  new_rot0( std::cos( pi/4 ) , std::sin( pi/4 ) ) ;
      qclab::QRotation< R >  new_rot1( std::cos( pi/6 ) , std::sin( pi/6 ) ) ;
      f3c::qgates::RotationXY< T >  Rxy2( new_rot0 , new_rot1 ) ;
      EXPECT_TRUE( Rxy == Rxy2 ) ;
      EXPECT_FALSE( Rxy != Rxy2 ) ;
      Rxy2.update( 1 , 2 ) ;
      EXPECT_TRUE( Rxy != Rxy2 ) ;
      EXPECT_FALSE( Rxy == Rxy2 ) ;
    }
  }


  //
  // constructors without qubits
  //
  {
    qclab::QRotation< R >         rot0( pi/4 ) ;
    qclab::QRotation< R >         rot1( pi/3 ) ;
    f3c::qgates::RotationXY< T >  Rxy( rot0 , rot1 ) ;

    EXPECT_EQ( Rxy.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rxy.fixed() ) ;              // fixed
    EXPECT_FALSE( Rxy.controlled() ) ;         // controlled

    auto qubits = Rxy.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;               // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;               // qubit1

    auto [ r0 , r1 ] = Rxy.rotations() ;
    EXPECT_TRUE( r0 == rot0 ) ;                // rot0
    EXPECT_TRUE( r1 == rot1 ) ;                // rot1

    auto [ theta0 , theta1 ] = Rxy.thetas() ;
    EXPECT_NEAR( theta0 , pi/4 , eps ) ;       // theta0
    EXPECT_NEAR( theta1 , pi/3 , eps ) ;       // theta1
  }

  {
    const R theta0( pi/4 ) ;
    const R theta1( pi/3 ) ;
    f3c::qgates::RotationXY< T >  Rxy( theta0 , theta1 ) ;

    EXPECT_EQ( Rxy.nbQubits() , 2 ) ;           // nbQubits
    EXPECT_FALSE( Rxy.fixed() ) ;               // fixed
    EXPECT_FALSE( Rxy.controlled() ) ;          // controlled

    auto qubits = Rxy.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;                // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;                // qubit1

    auto [ r0 , r1 ] = Rxy.rotations() ;
    EXPECT_NEAR( r0.theta() , theta0 , eps ) ;  // rot0
    EXPECT_NEAR( r1.theta() , theta1 , eps ) ;  // rot1

    auto [ th0 , th1 ] = Rxy.thetas() ;
    EXPECT_NEAR( th0 , theta0 , eps ) ;         // theta0
    EXPECT_NEAR( th1 , theta1 , eps ) ;         // theta1
  }


  //
  // constructors with qubit0 and qubit1
  //
  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    qclab::QRotation< R >         rot0( pi/4 ) ;
    qclab::QRotation< R >         rot1( pi/3 ) ;
    f3c::qgates::RotationXY< T >  Rxy( qubit0 , qubit1 , rot0 , rot1 ) ;

    EXPECT_EQ( Rxy.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rxy.fixed() ) ;              // fixed
    EXPECT_FALSE( Rxy.controlled() ) ;         // controlled

    auto qubits = Rxy.qubits() ;
    EXPECT_EQ( qubits[0] , qubit0 ) ;          // qubit0
    EXPECT_EQ( qubits[1] , qubit1 ) ;          // qubit1

    auto [ r0 , r1 ] = Rxy.rotations() ;
    EXPECT_TRUE( r0 == rot0 ) ;                // rot0
    EXPECT_TRUE( r1 == rot1 ) ;                // rot1

    auto [ theta0 , theta1 ] = Rxy.thetas() ;
    EXPECT_NEAR( theta0 , pi/4 , eps ) ;       // theta0
    EXPECT_NEAR( theta1 , pi/3 , eps ) ;       // theta1
  }

  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    const R theta0( pi/4 ) ;
    const R theta1( pi/3 ) ;
    f3c::qgates::RotationXY< T >  Rxy( qubit0 , qubit1 , theta0 , theta1 ) ;

    EXPECT_EQ( Rxy.nbQubits() , 2 ) ;           // nbQubits
    EXPECT_FALSE( Rxy.fixed() ) ;               // fixed
    EXPECT_FALSE( Rxy.controlled() ) ;          // controlled

    auto qubits = Rxy.qubits() ;
    EXPECT_EQ( qubits[0] , qubit0 ) ;           // qubit0
    EXPECT_EQ( qubits[1] , qubit1 ) ;           // qubit1

    auto [ r0 , r1 ] = Rxy.rotations() ;
    EXPECT_NEAR( r0.theta() , theta0 , eps ) ;  // rot0
    EXPECT_NEAR( r1.theta() , theta1 , eps ) ;  // rot1

    auto [ th0 , th1 ] = Rxy.thetas() ;
    EXPECT_NEAR( th0 , theta0 , eps ) ;         // theta0
    EXPECT_NEAR( th1 , theta1 , eps ) ;         // theta1
  }


  //
  // operators
  //
  {
    const R thetaA0 = pi/2 ;
    const R thetaA1 = pi/3 ;
    f3c::qgates::RotationXY< T >  RxyA( thetaA0 , thetaA1 ) ;

    const R thetaB0 = pi/4 ;
    const R thetaB1 = pi/5 ;
    f3c::qgates::RotationXY< T >  RxyB( thetaB0 , thetaB1 ) ;

    // operator *=
    RxyA *= RxyB ;
    const R theta0 = thetaA0 + thetaB0 ;
    const R theta1 = thetaA1 + thetaB1 ;
    auto [ th0 , th1 ] = RxyA.thetas() ;
    EXPECT_NEAR( th0 , theta0 , 10*eps ) ;    // thetaA0 = theta0
    EXPECT_NEAR( th1 , theta1 , 10*eps ) ;    // thetaA1 = theta1
    auto [ thB0 , thB1 ] = RxyB.thetas() ;
    EXPECT_NEAR( thB0 , thetaB0 , 10*eps ) ;  // thetaB0
    EXPECT_NEAR( thB1 , thetaB1 , 10*eps ) ;  // thetaB1
  }

  {
    const R thetaA0 = pi/2 ;
    const R thetaA1 = pi/3 ;
    f3c::qgates::RotationXY< T >  RxyA( thetaA0 , thetaA1 ) ;

    const R thetaB0 = pi/4 ;
    const R thetaB1 = pi/5 ;
    f3c::qgates::RotationXY< T >  RxyB( thetaB0 , thetaB1 ) ;

    // operator /=
    RxyA /= RxyB ;
    const R theta0 = thetaA0 - thetaB0 ;
    const R theta1 = thetaA1 - thetaB1 ;
    auto [ th0 , th1 ] = RxyA.thetas() ;
    EXPECT_NEAR( th0 , theta0 , 10*eps ) ;    // thetaA0 = theta0
    EXPECT_NEAR( th1 , theta1 , 10*eps ) ;    // thetaA1 = theta1
    auto [ thB0 , thB1 ] = RxyB.thetas() ;
    EXPECT_NEAR( thB0 , thetaB0 , 10*eps ) ;  // thetaB0
    EXPECT_NEAR( thB1 , thetaB1 , 10*eps ) ;  // thetaB1
  }

  {
    const R thetaA0 = pi/2 ;
    const R thetaA1 = pi/3 ;
    f3c::qgates::RotationXY< T >  RxyA( thetaA0 , thetaA1 ) ;

    const R thetaB0 = pi/4 ;
    const R thetaB1 = pi/5 ;
    f3c::qgates::RotationXY< T >  RxyB( thetaB0 , thetaB1 ) ;

    // operator *
    f3c::qgates::RotationXY< T > RxyC = RxyA * RxyB ;
    const R thetaC0 = thetaA0 + thetaB0 ;
    const R thetaC1 = thetaA1 + thetaB1 ;
    auto [ thA0 , thA1 ] = RxyA.thetas() ;
    EXPECT_NEAR( thA0 , thetaA0 , 10*eps ) ;  // thetaA0
    EXPECT_NEAR( thA1 , thetaA1 , 10*eps ) ;  // thetaA1
    auto [ thB0 , thB1 ] = RxyB.thetas() ;
    EXPECT_NEAR( thB0 , thetaB0 , 10*eps ) ;  // thetaB0
    EXPECT_NEAR( thB1 , thetaB1 , 10*eps ) ;  // thetaB1
    auto [ thC0 , thC1 ] = RxyC.thetas() ;
    EXPECT_NEAR( thC0 , thetaC0 , 10*eps ) ;  // thetaC0
    EXPECT_NEAR( thC1 , thetaC1 , 10*eps ) ;  // thetaC1
  }

  {
    const R thetaA0 = pi/2 ;
    const R thetaA1 = pi/3 ;
    f3c::qgates::RotationXY< T >  RxyA( thetaA0 , thetaA1 ) ;

    const R thetaB0 = pi/4 ;
    const R thetaB1 = pi/5 ;
    f3c::qgates::RotationXY< T >  RxyB( thetaB0 , thetaB1 ) ;

    // operator /
    f3c::qgates::RotationXY< T > RxyC = RxyA / RxyB ;
    const R thetaC0 = thetaA0 - thetaB0 ;
    const R thetaC1 = thetaA1 - thetaB1 ;
    auto [ thA0 , thA1 ] = RxyA.thetas() ;
    EXPECT_NEAR( thA0 , thetaA0 , 10*eps ) ;  // thetaA0
    EXPECT_NEAR( thA1 , thetaA1 , 10*eps ) ;  // thetaA1
    auto [ thB0 , thB1 ] = RxyB.thetas() ;
    EXPECT_NEAR( thB0 , thetaB0 , 10*eps ) ;  // thetaB0
    EXPECT_NEAR( thB1 , thetaB1 , 10*eps ) ;  // thetaB1
    auto [ thC0 , thC1 ] = RxyC.thetas() ;
    EXPECT_NEAR( thC0 , thetaC0 , 10*eps ) ;  // thetaC0
    EXPECT_NEAR( thC1 , thetaC1 , 10*eps ) ;  // thetaC1
  }

  {
    const R thetaA0 = pi/2 ;
    const R thetaA1 = pi/3 ;
    f3c::qgates::RotationXY< T >  RxyA( thetaA0 , thetaA1 ) ;

    // inv
    f3c::qgates::RotationXY< T >  RxyB = RxyA.inv() ;
    auto [ thA0 , thA1 ] = RxyA.thetas() ;
    EXPECT_NEAR( thA0 , thetaA0 , 10*eps ) ;  // thetaA0
    EXPECT_NEAR( thA1 , thetaA1 , 10*eps ) ;  // thetaA1
    auto [ thB0 , thB1 ] = RxyB.thetas() ;
    EXPECT_NEAR( thB0 , -thetaA0 , 10*eps ) ;  // thetaB0
    EXPECT_NEAR( thB1 , -thetaA1 , 10*eps ) ;  // thetaB1
  }

}


/*
 * complex float
 */
TEST( f3c_qgates_RotationXY , complex_float ) {
  test_f3c_qgates_RotationXY< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( f3c_qgates_RotationXY , complex_double ) {
  test_f3c_qgates_RotationXY< std::complex< double > >() ;
}

