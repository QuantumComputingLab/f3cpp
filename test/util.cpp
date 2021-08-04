#include <gtest/gtest.h>
#include "f3c/util.hpp"

template <typename R>
void test_f3c_util_eig22() {

  using T = std::complex< R > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  std::array< T , 4 >  A = { 2 , 0 , 0 , 3 } ;
  std::array< T , 4 >  B = { 1 , 0 , 0 , 1 } ;

  {
    T a1 ;
    T a2 ;
    T b ;
    f3c::eig22( A.data() , B.data() , a1 , a2 , b ) ;

    EXPECT_EQ( a1 , T(3) ) ;
    EXPECT_EQ( a2 , T(2) ) ;
    EXPECT_EQ( b  , T(1) ) ;
  }

  A[0] = T( -6.165451364571467e-01 , -3.168692408037940e-01 ) ;
  A[1] = T( -8.745993399503420e-02 ,  7.154136837446220e-01 ) ;
  A[2] = T( -3.345998130349203e-01 ,  6.383642722305675e-01 ) ;
  A[3] = T( -6.873841671402685e-01 , -8.964947195278750e-02 ) ;

  B[0] = T( -7.335269870327210e-01 ,  3.507737419764394e-01 ) ;
  B[1] = T(  3.839162873894437e-01 , -4.376119576881443e-01 ) ;
  B[2] = T( -5.670896578574434e-01 , -1.315494628862925e-01 ) ;
  B[3] = T( -8.004273135165704e-01 ,  1.428991761421370e-01 ) ;

  {
    T a1 ;
    T a2 ;
    T b ;
    f3c::eig22( A.data() , B.data() , a1 , a2 , b ) ;

    EXPECT_NEAR( std::real(a1) ,  3.983405377897217e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag(a1) ,  9.172376005994276e-01 , 10*eps ) ;
    EXPECT_NEAR( std::real(a2) ,  3.814934699105096e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag(a2) , -9.243715337544955e-01 , 10*eps ) ;
    EXPECT_NEAR( std::real(b)  ,  8.122923309298512e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag(b)  , -5.832505200276707e-01 , 10*eps ) ;
  }

}


template <typename R>
void test_f3c_util_rotateToZero() {

  using T = std::complex< R > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    T x = 1 ;
    T y = 0 ;
    auto [ c , s ] = f3c::rotateToZero( x , y ) ;

    EXPECT_NEAR( std::real(c) , 1 , 10*eps ) ;
    EXPECT_NEAR( std::imag(c) , 0 , 10*eps ) ;
    EXPECT_NEAR( std::real(s) , 0 , 10*eps ) ;
    EXPECT_NEAR( std::imag(s) , 0 , 10*eps ) ;
  }

  {
    T x = 0 ;
    T y = 1 ;
    auto [ c , s ] = f3c::rotateToZero( x , y ) ;

    EXPECT_NEAR( std::real(c) , 0 , 10*eps ) ;
    EXPECT_NEAR( std::imag(c) , 0 , 10*eps ) ;
    EXPECT_NEAR( std::real(s) , 1 , 10*eps ) ;
    EXPECT_NEAR( std::imag(s) , 0 , 10*eps ) ;
  }

  {
    T x( 3.277245845061543e-01 , -8.879587326190993e-01 ) ;
    T y( 3.212070455035705e-01 , -3.085319747803702e-02 ) ;
    auto [ c , s ] = f3c::rotateToZero( x , y ) ;

    EXPECT_NEAR( std::real(c) , 3.277245845061542e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag(c) , 8.879587326190992e-01 , 10*eps ) ;
    EXPECT_NEAR( std::real(s) , 3.212070455035705e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag(s) , 3.085319747803701e-02 , 10*eps ) ;
  }

}


template <typename R>
void test_f3c_util_rotateToZeroL() {

  using T = std::complex< R > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    T x( 3.277245845061543e-01 , -8.879587326190993e-01 ) ;
    T y( 3.212070455035705e-01 , -3.085319747803702e-02 ) ;
    std::array< T , 4 >  G ;
    f3c::rotateToZeroL( x , y , G.data() ) ;

    EXPECT_NEAR( std::real( G[0] ) ,  3.277245845061542e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[0] ) ,  8.879587326190992e-01 , 10*eps ) ;
    EXPECT_NEAR( std::real( G[1] ) , -3.212070455035705e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[1] ) ,  3.085319747803701e-02 , 10*eps ) ;
    EXPECT_NEAR( std::real( G[2] ) ,  3.212070455035705e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[2] ) ,  3.085319747803701e-02 , 10*eps ) ;
    EXPECT_NEAR( std::real( G[3] ) ,  3.277245845061542e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[3] ) , -8.879587326190992e-01 , 10*eps ) ;
  }

}


template <typename R>
void test_f3c_util_rotateToZeroR() {

  using T = std::complex< R > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    T x( 3.277245845061543e-01 , -8.879587326190993e-01 ) ;
    T y( 3.212070455035705e-01 , -3.085319747803702e-02 ) ;
    std::array< T , 4 >  G ;
    f3c::rotateToZeroR( x , y , G.data() ) ;

    EXPECT_NEAR( std::real( G[0] ) ,  3.277245845061542e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[0] ) , -8.879587326190992e-01 , 10*eps ) ;
    EXPECT_NEAR( std::real( G[1] ) , -3.212070455035705e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[1] ) ,  3.085319747803701e-02 , 10*eps ) ;
    EXPECT_NEAR( std::real( G[2] ) ,  3.212070455035705e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[2] ) ,  3.085319747803701e-02 , 10*eps ) ;
    EXPECT_NEAR( std::real( G[3] ) ,  3.277245845061542e-01 , 10*eps ) ;
    EXPECT_NEAR( std::imag( G[3] ) ,  8.879587326190992e-01 , 10*eps ) ;
  }

}


template <typename R>
void test_f3c_util_diagonalize22() {

  using T = std::complex< R > ;
  using M = std::array< T , 4 > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  M  A( { T( -1.773751566188252e-01 ,  1.978110534643607e-01 ) ,
          T( -1.960534878073328e-01 ,  1.587699089974059e+00 ) ,
          T(  1.419310150642549e+00 , -8.044659563495471e-01 ) ,
          T(  2.915843739841825e-01 ,  6.966244158496073e-01 ) } ) ;
  M  B( { T(  2.915843739841825e-01 , -6.966244158496073e-01 ) ,
          T( -1.419310150642549e+00 , -8.044659563495471e-01 ) ,
          T(  1.960534878073328e-01 ,  1.587699089974059e+00 ) ,
          T( -1.773751566188252e-01 , -1.978110534643607e-01 ) } ) ;
  M  Q ;
  M  Z ;

  f3c::diagonalize22( A.data() , B.data() , Q.data() , Z.data() ) ;

  EXPECT_NEAR( std::real( A[0] ) , -5.978906623049421e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( A[0] ) , -1.246098865369929e+00 , 10*eps ) ;
  EXPECT_EQ( A[1] , T(0) ) ;
  EXPECT_EQ( A[2] , T(0) ) ;
  EXPECT_NEAR( std::real( A[3] ) ,  1.987836689883845e+00 , 10*eps ) ;
  EXPECT_NEAR( std::imag( A[3] ) ,  0.000000000000000e+00 , 10*eps ) ;

  EXPECT_NEAR( std::real( B[0] ) ,  1.987836689883845e+00 , 10*eps ) ;
  EXPECT_NEAR( std::imag( B[0] ) , -5.551115123125783e-17 , 10*eps ) ;
  EXPECT_EQ( B[1] , T(0) ) ;
  EXPECT_EQ( B[2] , T(0) ) ;
  EXPECT_NEAR( std::real( B[3] ) , -5.978906623049423e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( B[3] ) ,  1.246098865369929e+00 , 10*eps ) ;

  EXPECT_NEAR( std::real( Q[0] ) , -4.383880931972775e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Q[0] ) , -6.365074287717316e-01 , 10*eps ) ;
  EXPECT_NEAR( std::real( Q[1] ) , -5.845545419235243e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Q[1] ) , -2.469213647658552e-01 , 10*eps ) ;
  EXPECT_NEAR( std::real( Q[2] ) ,  5.845545419235243e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Q[2] ) , -2.469213647658552e-01 , 10*eps ) ;
  EXPECT_NEAR( std::real( Q[3] ) , -4.383880931972775e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Q[3] ) ,  6.365074287717316e-01 , 10*eps ) ;

  EXPECT_NEAR( std::real( Z[0] ) , -8.046625602671207e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Z[0] ) ,  4.188811247454753e-17 , 10*eps ) ;
  EXPECT_NEAR( std::real( Z[1] ) ,  3.884149668603150e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Z[1] ) ,  4.490567643664469e-01 , 10*eps ) ;
  EXPECT_NEAR( std::real( Z[2] ) , -3.884149668603150e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Z[2] ) ,  4.490567643664469e-01 , 10*eps ) ;
  EXPECT_NEAR( std::real( Z[3] ) , -8.046625602671207e-01 , 10*eps ) ;
  EXPECT_NEAR( std::imag( Z[3] ) , -4.188811247454753e-17 , 10*eps ) ;

}


template <typename R>
void test_f3c_util() {

  test_f3c_util_eig22< R >() ;
  test_f3c_util_rotateToZero< R >() ;
  test_f3c_util_rotateToZeroL< R >() ;
  test_f3c_util_rotateToZeroR< R >() ;
  test_f3c_util_diagonalize22< R >() ;

}


/*
 * float
 */
TEST( f3c_util , float ) {
  test_f3c_util< float >() ;
}

/*
 * double
 */
TEST( f3c_util , double ) {
  test_f3c_util< double >() ;
}

