#include <gtest/gtest.h>
#include "f3c/qgates/RotationTFXYMatrix.hpp"
#include "f3c/qgates/RotationTFXY.hpp"

template <typename T>
void test_f3c_qgates_RotationTFXYMatrix() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    f3c::qgates::RotationTFXYMatrix< T >  TFRxy ;

    EXPECT_EQ( TFRxy.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_FALSE( TFRxy.fixed() ) ;       // fixed
    EXPECT_FALSE( TFRxy.controlled() ) ;  // controlled

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( TFRxy.matrix() == eye ) ;

    // qubits
    EXPECT_EQ( TFRxy.qubits().size() , 2 ) ;
    EXPECT_EQ( TFRxy.qubits()[0] , 0 ) ;
    EXPECT_EQ( TFRxy.qubits()[1] , 1 ) ;
    int qnew[2] = { 3 , 5 } ;
    TFRxy.setQubits( &qnew[0] ) ;
    EXPECT_EQ( TFRxy.qubits().size() , 2 ) ;
    EXPECT_EQ( TFRxy.qubits()[0] , 3 ) ;
    EXPECT_EQ( TFRxy.qubits()[1] , 5 ) ;

    // update(a,b,c,d)
    const T a( std::sqrt( 2. ) / 2 ) ;
    const T b( 0 , std::sqrt( 2. ) / 2 ) ;
    const T c( -std::sqrt( 2. ) / 2 ) ;
    const T d( 0 , -std::sqrt( 2. ) / 2 ) ;
    TFRxy.update( a , b , c , d ) ;
    EXPECT_EQ( TFRxy.a() , a ) ;
    EXPECT_EQ( TFRxy.b() , b ) ;
    EXPECT_EQ( TFRxy.c() , c ) ;
    EXPECT_EQ( TFRxy.d() , d ) ;

    // matrix
    EXPECT_EQ( TFRxy.matrix()(0,0) ,   a  ) ;
    EXPECT_EQ( TFRxy.matrix()(1,0) , T(0) ) ;
    EXPECT_EQ( TFRxy.matrix()(2,0) , T(0) ) ;
    EXPECT_EQ( TFRxy.matrix()(3,0) ,   d  ) ;

    EXPECT_EQ( TFRxy.matrix()(0,1) , T(0) ) ;
    EXPECT_EQ( TFRxy.matrix()(1,1) ,   b  ) ;
    EXPECT_EQ( TFRxy.matrix()(2,1) ,   c  ) ;
    EXPECT_EQ( TFRxy.matrix()(3,1) , T(0) ) ;

    EXPECT_EQ( TFRxy.matrix()(0,2) , T(0) ) ;
    EXPECT_EQ( TFRxy.matrix()(1,2) , -std::conj( c ) ) ;
    EXPECT_EQ( TFRxy.matrix()(2,2) ,  std::conj( b ) ) ;
    EXPECT_EQ( TFRxy.matrix()(3,2) , T(0) ) ;

    EXPECT_EQ( TFRxy.matrix()(0,3) , -std::conj( d ) ) ;
    EXPECT_EQ( TFRxy.matrix()(1,3) , T(0) ) ;
    EXPECT_EQ( TFRxy.matrix()(2,3) , T(0) ) ;
    EXPECT_EQ( TFRxy.matrix()(3,3) ,  std::conj( a ) ) ;

    // print
    TFRxy.print() ;

    // toQASM
    {
      std::stringstream qasm ;
      EXPECT_EQ( TFRxy.toQASM( qasm ) , 0 ) ;
      std::stringstream qasm_check ;
      f3c::qgates::RotationTFXY< T >  tmp( TFRxy ) ;
      tmp.toQASM( qasm_check ) ;
      EXPECT_EQ( qasm.str() , qasm_check.str() ) ;
    }

    // operators == and !=
    {
      f3c::qgates::RotationTFXYMatrix< T >  TFRxy2( a , b , c , d ) ;
      EXPECT_TRUE( TFRxy == TFRxy2 ) ;
      EXPECT_FALSE( TFRxy != TFRxy2 ) ;
      TFRxy2.update( 1 , 1 , 0 , 0 ) ;
      EXPECT_TRUE( TFRxy != TFRxy2 ) ;
      EXPECT_FALSE( TFRxy == TFRxy2 ) ;
    }
  }


  //
  // constructor without qubits
  //
  {
    const T a( std::sqrt( 2. ) / 2 ) ;
    const T b( 0 , std::sqrt( 2. ) / 2 ) ;
    const T c( -std::sqrt( 2. ) / 2 ) ;
    const T d( 0 , -std::sqrt( 2. ) / 2 ) ;
    f3c::qgates::RotationTFXYMatrix< T >  TFRxy( a , b , c , d ) ;

    EXPECT_EQ( TFRxy.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_FALSE( TFRxy.fixed() ) ;       // fixed
    EXPECT_FALSE( TFRxy.controlled() ) ;  // controlled

    auto qubits = TFRxy.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;          // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;          // qubit1

    EXPECT_EQ( TFRxy.a() , a ) ;          // a
    EXPECT_EQ( TFRxy.b() , b ) ;          // b
    EXPECT_EQ( TFRxy.c() , c ) ;          // c
    EXPECT_EQ( TFRxy.d() , d ) ;          // d
  }


  //
  // constructor with qubit0 and qubit1
  //
  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    const T a( std::sqrt( 2. ) / 2 ) ;
    const T b( 0 , std::sqrt( 2. ) / 2 ) ;
    const T c( -std::sqrt( 2. ) / 2 ) ;
    const T d( 0 , -std::sqrt( 2. ) / 2 ) ;
    f3c::qgates::RotationTFXYMatrix< T >  TFRxy( qubit0 , qubit1 ,
                                                 a , b , c , d ) ;

    EXPECT_EQ( TFRxy.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_FALSE( TFRxy.fixed() ) ;       // fixed
    EXPECT_FALSE( TFRxy.controlled() ) ;  // controlled

    auto qubits = TFRxy.qubits() ;
    EXPECT_EQ( qubits[0] , qubit0 ) ;     // qubit0
    EXPECT_EQ( qubits[1] , qubit1 ) ;     // qubit1

    EXPECT_EQ( TFRxy.a() , a ) ;          // a
    EXPECT_EQ( TFRxy.b() , b ) ;          // b
    EXPECT_EQ( TFRxy.c() , c ) ;          // c
    EXPECT_EQ( TFRxy.d() , d ) ;          // d
  }


  //
  // constructor with gate
  //
  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    const R theta0 =  pi/3 ;
    const R theta1 = -0.33 ;
    const R theta2 =  pi/2 ;
    const R theta3 = -pi/7 ;
    const R theta4 = -pi/5 ;
    const R theta5 =  pi/9 ;
    f3c::qgates::RotationTFXY< T >  TFRxy( qubit0 , qubit1 ,
                                           theta0 , theta1 , theta2 ,
                                           theta3 , theta4 , theta5 ) ;

    f3c::qgates::RotationTFXYMatrix< T >  TFRxyMat( TFRxy ) ;

    EXPECT_EQ( TFRxyMat.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_FALSE( TFRxyMat.fixed() ) ;       // fixed
    EXPECT_FALSE( TFRxyMat.controlled() ) ;  // controlled

    auto qubits = TFRxyMat.qubits() ;
    EXPECT_EQ( qubits[0] , qubit0 ) ;        // qubit0
    EXPECT_EQ( qubits[1] , qubit1 ) ;        // qubit1

    EXPECT_EQ( TFRxyMat.a() , TFRxy.a() ) ;  // a
    EXPECT_EQ( TFRxyMat.b() , TFRxy.b() ) ;  // b
    EXPECT_EQ( TFRxyMat.c() , TFRxy.c() ) ;  // c
    EXPECT_EQ( TFRxyMat.d() , TFRxy.d() ) ;  // d
  }


  //
  // operators
  //
  {
    const T a1 = T( -5.049262824994509e-01 , -1.183280215457461e-01 ) ;
    const T b1 = T( -1.646607339353139e-01 , -8.390082069075593e-01 ) ;
    const T c1 = T( -2.653946050859590e-01 ,  4.455533357892336e-01 ) ;
    const T d1 = T( -4.895245818757411e-01 , -7.010089958750423e-01 ) ;

    const T a2 = T(  3.672139606908904e-01 ,  5.490118451599342e-01 ) ;
    const T b2 = T( -6.062117964174196e-02 ,  8.916968089208992e-01 ) ;
    const T c2 = T( -3.703676033540890e-01 ,  2.530409293471655e-01 ) ;
    const T d2 = T(  4.740375143541831e-02 , -7.493282226752457e-01 ) ;

    f3c::qgates::RotationTFXYMatrix< T >  RxyA( a1 , b1 , c1 , d1 ) ;
    f3c::qgates::RotationTFXYMatrix< T >  RxyB( a2 , b2 , c2 , d2 ) ;

    // operator *=
    RxyA *= RxyB ;

    EXPECT_NEAR( std::real( RxyA.a() ) , -6.225330179682356e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.a() ) ,  7.938282960994436e-02 , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyA.b() ) ,  5.470860746991210e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.b() ) ,  1.897039855149413e-03 , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyA.c() ) ,  6.866754395581218e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.c() ) ,  4.787171072959016e-01 , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyA.d() ) , -6.772244289315784e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.d() ) ,  3.840808258176608e-01 , 10*eps ) ;

    EXPECT_NEAR( std::real( RxyB.a() ) , std::real( a2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.a() ) , std::imag( a2 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyB.b() ) , std::real( b2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.b() ) , std::imag( b2 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyB.c() ) , std::real( c2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.c() ) , std::imag( c2 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyB.d() ) , std::real( d2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.d() ) , std::imag( d2 ) , 10*eps ) ;
  }

  {
    const T a1 = T( -5.049262824994509e-01 , -1.183280215457461e-01 ) ;
    const T b1 = T( -1.646607339353139e-01 , -8.390082069075593e-01 ) ;
    const T c1 = T( -2.653946050859590e-01 ,  4.455533357892336e-01 ) ;
    const T d1 = T( -4.895245818757411e-01 , -7.010089958750423e-01 ) ;

    const T a2 = T(  3.672139606908904e-01 ,  5.490118451599342e-01 ) ;
    const T b2 = T( -6.062117964174196e-02 ,  8.916968089208992e-01 ) ;
    const T c2 = T( -3.703676033540890e-01 ,  2.530409293471655e-01 ) ;
    const T d2 = T(  4.740375143541831e-02 , -7.493282226752457e-01 ) ;

    f3c::qgates::RotationTFXYMatrix< T >  RxyA( a1 , b1 , c1 , d1 ) ;
    f3c::qgates::RotationTFXYMatrix< T >  RxyB( a2 , b2 , c2 , d2 ) ;

    // operator *
    f3c::qgates::RotationTFXYMatrix< T > RxyC = RxyA * RxyB ;

    EXPECT_NEAR( std::real( RxyA.a() ) , std::real( a1 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.a() ) , std::imag( a1 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyA.b() ) , std::real( b1 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.b() ) , std::imag( b1 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyA.c() ) , std::real( c1 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.c() ) , std::imag( c1 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyA.d() ) , std::real( d1 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyA.d() ) , std::imag( d1 ) , 10*eps ) ;

    EXPECT_NEAR( std::real( RxyB.a() ) , std::real( a2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.a() ) , std::imag( a2 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyB.b() ) , std::real( b2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.b() ) , std::imag( b2 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyB.c() ) , std::real( c2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.c() ) , std::imag( c2 ) , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyB.d() ) , std::real( d2 ) , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyB.d() ) , std::imag( d2 ) , 10*eps ) ;

    EXPECT_NEAR( std::real( RxyC.a() ) , -6.225330179682356e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyC.a() ) ,  7.938282960994436e-02 , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyC.b() ) ,  5.470860746991210e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyC.b() ) ,  1.897039855149413e-03 , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyC.c() ) ,  6.866754395581218e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyC.c() ) ,  4.787171072959016e-01 , 10*eps ) ;
    EXPECT_NEAR( std::real( RxyC.d() ) , -6.772244289315784e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( RxyC.d() ) ,  3.840808258176608e-01 , 10*eps ) ;
  }

}


/*
 * complex float
 */
TEST( f3c_qgates_RotationTFXYMatrix , complex_float ) {
  test_f3c_qgates_RotationTFXYMatrix< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( f3c_qgates_RotationTFXYMatrix , complex_double ) {
  test_f3c_qgates_RotationTFXYMatrix< std::complex< double > >() ;
}

