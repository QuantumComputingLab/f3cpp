#include <gtest/gtest.h>
#include "f3c/concepts.hpp"
#include "f3c/qgates/RotationTFXYMatrix.hpp"

template <typename T>
void test_f3c_concepts_is_one_axis() {

  using X = qclab::qgates::RotationX< T > ;
  using Y = qclab::qgates::RotationY< T > ;
  using Z = qclab::qgates::RotationZ< T > ;

  using XX = qclab::qgates::RotationXX< T > ;
  using YY = qclab::qgates::RotationYY< T > ;
  using ZZ = qclab::qgates::RotationZZ< T > ;

  using XY = f3c::qgates::RotationXY< T > ;
  using XZ = f3c::qgates::RotationXZ< T > ;
  using YZ = f3c::qgates::RotationYZ< T > ;

  EXPECT_FALSE( f3c::is_one_axis1< T >::value ) ;
  EXPECT_FALSE( f3c::is_one_axis1_v< T > ) ;

  EXPECT_TRUE( f3c::is_one_axis1_v< X > ) ;
  EXPECT_TRUE( f3c::is_one_axis1_v< Y > ) ;
  EXPECT_TRUE( f3c::is_one_axis1_v< Z > ) ;

  EXPECT_FALSE( f3c::is_one_axis1_v< XX > ) ;
  EXPECT_FALSE( f3c::is_one_axis1_v< YY > ) ;
  EXPECT_FALSE( f3c::is_one_axis1_v< ZZ > ) ;

  EXPECT_FALSE( f3c::is_one_axis1_v< XY > ) ;
  EXPECT_FALSE( f3c::is_one_axis1_v< XZ > ) ;
  EXPECT_FALSE( f3c::is_one_axis1_v< YZ > ) ;


  EXPECT_FALSE( f3c::is_one_axis2< T >::value ) ;
  EXPECT_FALSE( f3c::is_one_axis2_v< T > ) ;

  EXPECT_FALSE( f3c::is_one_axis2_v< X > ) ;
  EXPECT_FALSE( f3c::is_one_axis2_v< Y > ) ;
  EXPECT_FALSE( f3c::is_one_axis2_v< Z > ) ;

  EXPECT_TRUE( f3c::is_one_axis2_v< XX > ) ;
  EXPECT_TRUE( f3c::is_one_axis2_v< YY > ) ;
  EXPECT_TRUE( f3c::is_one_axis2_v< ZZ > ) ;

  EXPECT_FALSE( f3c::is_one_axis2_v< XY > ) ;
  EXPECT_FALSE( f3c::is_one_axis2_v< XZ > ) ;
  EXPECT_FALSE( f3c::is_one_axis2_v< YZ > ) ;


  EXPECT_FALSE( f3c::is_one_axis_v< T > ) ;

  EXPECT_TRUE( f3c::is_one_axis_v< X > ) ;
  EXPECT_TRUE( f3c::is_one_axis_v< Y > ) ;
  EXPECT_TRUE( f3c::is_one_axis_v< Z > ) ;

  EXPECT_TRUE( f3c::is_one_axis_v< XX > ) ;
  EXPECT_TRUE( f3c::is_one_axis_v< YY > ) ;
  EXPECT_TRUE( f3c::is_one_axis_v< ZZ > ) ;

  EXPECT_FALSE( f3c::is_one_axis_v< XY > ) ;
  EXPECT_FALSE( f3c::is_one_axis_v< XZ > ) ;
  EXPECT_FALSE( f3c::is_one_axis_v< YZ > ) ;

}

template <typename T>
void test_f3c_concepts_is_same_axis() {

  using X = qclab::qgates::RotationX< T > ;
  using Y = qclab::qgates::RotationY< T > ;
  using Z = qclab::qgates::RotationZ< T > ;

  using XX = qclab::qgates::RotationXX< T > ;
  using YY = qclab::qgates::RotationYY< T > ;
  using ZZ = qclab::qgates::RotationZZ< T > ;

  EXPECT_FALSE( ( f3c::is_same_axis< T , T >::value ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< X , X > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< Y , Y > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< Z , Z > ) ) ;

  EXPECT_FALSE( ( f3c::is_same_axis_v< X , Y > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< Y , X > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< X , Z > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< Z , X > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< Y , Z > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< Z , Y > ) ) ;

  EXPECT_TRUE( ( f3c::is_same_axis_v< XX , XX > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< YY , YY > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< ZZ , ZZ > ) ) ;

  EXPECT_FALSE( ( f3c::is_same_axis_v< XX , YY > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< YY , XX > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< XX , ZZ > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< ZZ , XX > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< YY , ZZ > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< ZZ , YY > ) ) ;

  EXPECT_TRUE( ( f3c::is_same_axis_v< X , XX > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< XX , X > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< Y , YY > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< YY , Y > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< Z , ZZ > ) ) ;
  EXPECT_TRUE( ( f3c::is_same_axis_v< ZZ , Z > ) ) ;

  EXPECT_FALSE( ( f3c::is_same_axis_v< X , YY > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< YY , X > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< Y , XX > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< XX , Y > ) ) ;

  EXPECT_FALSE( ( f3c::is_same_axis_v< X , ZZ > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< ZZ , X > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< Z , XX > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< XX , Z > ) ) ;

  EXPECT_FALSE( ( f3c::is_same_axis_v< Y , ZZ > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< ZZ , Y > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< Z , YY > ) ) ;
  EXPECT_FALSE( ( f3c::is_same_axis_v< YY , Z > ) ) ;

}

template <typename T>
void test_f3c_concepts_is_two_axes() {

  using X = qclab::qgates::RotationX< T > ;
  using Y = qclab::qgates::RotationY< T > ;
  using Z = qclab::qgates::RotationZ< T > ;

  using XX = qclab::qgates::RotationXX< T > ;
  using YY = qclab::qgates::RotationYY< T > ;
  using ZZ = qclab::qgates::RotationZZ< T > ;

  using XY = f3c::qgates::RotationXY< T > ;
  using XZ = f3c::qgates::RotationXZ< T > ;
  using YZ = f3c::qgates::RotationYZ< T > ;

  EXPECT_FALSE( f3c::is_two_axes< T >::value ) ;
  EXPECT_FALSE( f3c::is_two_axes_v< T > ) ;

  EXPECT_FALSE( f3c::is_two_axes_v< X > ) ;
  EXPECT_FALSE( f3c::is_two_axes_v< Y > ) ;
  EXPECT_FALSE( f3c::is_two_axes_v< Z > ) ;

  EXPECT_FALSE( f3c::is_two_axes_v< XX > ) ;
  EXPECT_FALSE( f3c::is_two_axes_v< YY > ) ;
  EXPECT_FALSE( f3c::is_two_axes_v< ZZ > ) ;

  EXPECT_TRUE( f3c::is_two_axes_v< XY > ) ;
  EXPECT_TRUE( f3c::is_two_axes_v< XZ > ) ;
  EXPECT_TRUE( f3c::is_two_axes_v< YZ > ) ;

}

template <typename T>
void test_f3c_concepts_is_TF_two_axes() {

  using X = qclab::qgates::RotationX< T > ;
  using Y = qclab::qgates::RotationY< T > ;
  using Z = qclab::qgates::RotationZ< T > ;

  using XX = qclab::qgates::RotationXX< T > ;
  using YY = qclab::qgates::RotationYY< T > ;
  using ZZ = qclab::qgates::RotationZZ< T > ;

  using XY = f3c::qgates::RotationXY< T > ;
  using XZ = f3c::qgates::RotationXZ< T > ;
  using YZ = f3c::qgates::RotationYZ< T > ;

  using TFXY = f3c::qgates::RotationTFXY< T > ;
  using TFXZ = f3c::qgates::RotationTFXZ< T > ;
  using TFYZ = f3c::qgates::RotationTFYZ< T > ;

  using TFXYM = f3c::qgates::RotationTFXYMatrix< T > ;

  EXPECT_FALSE( f3c::is_TF_two_axes< T >::value ) ;
  EXPECT_FALSE( f3c::is_TF_two_axes_v< T > ) ;

  EXPECT_FALSE( f3c::is_TF_two_axes_v< X > ) ;
  EXPECT_FALSE( f3c::is_TF_two_axes_v< Y > ) ;
  EXPECT_FALSE( f3c::is_TF_two_axes_v< Z > ) ;

  EXPECT_FALSE( f3c::is_TF_two_axes_v< XX > ) ;
  EXPECT_FALSE( f3c::is_TF_two_axes_v< YY > ) ;
  EXPECT_FALSE( f3c::is_TF_two_axes_v< ZZ > ) ;

  EXPECT_FALSE( f3c::is_TF_two_axes_v< XY > ) ;
  EXPECT_FALSE( f3c::is_TF_two_axes_v< XZ > ) ;
  EXPECT_FALSE( f3c::is_TF_two_axes_v< YZ > ) ;

  EXPECT_TRUE( f3c::is_TF_two_axes_v< TFXY > ) ;
  EXPECT_TRUE( f3c::is_TF_two_axes_v< TFXZ > ) ;
  EXPECT_TRUE( f3c::is_TF_two_axes_v< TFYZ > ) ;

  EXPECT_FALSE( f3c::is_TF_two_axes_v< TFXYM > ) ;

}


template <typename R>
void test_f3c_concepts() {

  test_f3c_concepts_is_one_axis< R >() ;
  test_f3c_concepts_is_same_axis< R >() ;
  test_f3c_concepts_is_two_axes< R >() ;
  test_f3c_concepts_is_TF_two_axes< R >() ;

}


/*
 * complex float
 */
TEST( f3c_concepts , complex_float ) {
  test_f3c_concepts< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( f3c_concepts , complex_double ) {
  test_f3c_concepts< std::complex< double > >() ;
}

