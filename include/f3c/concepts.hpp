//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_concepts_hpp
#define f3c_concepts_hpp

#include "qclab/qgates/RotationX.hpp"
#include "qclab/qgates/RotationY.hpp"
#include "qclab/qgates/RotationZ.hpp"
#include "qclab/qgates/RotationXX.hpp"
#include "qclab/qgates/RotationYY.hpp"
#include "qclab/qgates/RotationZZ.hpp"
#include "f3c/qgates/RotationXY.hpp"
#include "f3c/qgates/RotationXZ.hpp"
#include "f3c/qgates/RotationYZ.hpp"
#include "f3c/qgates/RotationTFXY.hpp"
#include "f3c/qgates/RotationTFXZ.hpp"
#include "f3c/qgates/RotationTFYZ.hpp"
#include <type_traits>
//#include <concepts> // TODO: c++20

namespace f3c {

  /// Default helper class for is_one_axis1.
  template <typename T>
  struct is_one_axis1_helper
  : std::false_type { } ;

  /// Template specialized helper class for is_one_axis1.
  template <typename T>
  struct is_one_axis1_helper< qclab::qgates::RotationX< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_one_axis1.
  template <typename T>
  struct is_one_axis1_helper< qclab::qgates::RotationY< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_one_axis1.
  template <typename T>
  struct is_one_axis1_helper< qclab::qgates::RotationZ< T > >
  : std::true_type { } ;

  /// Returns type of is_one_axis1.
  template <typename T>
  struct is_one_axis1
  : is_one_axis1_helper< T >::type { } ;

  /// Checks if T is a one axis rotation gate acting on one qubit.
  template <typename T>
  inline constexpr bool is_one_axis1_v = is_one_axis1< T >::value ;


  /// Default helper class for is_one_axis2.
  template <typename T>
  struct is_one_axis2_helper
  : std::false_type { } ;

  /// Template specialized helper class for is_one_axis2.
  template <typename T>
  struct is_one_axis2_helper< qclab::qgates::RotationXX< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_one_axis2.
  template <typename T>
  struct is_one_axis2_helper< qclab::qgates::RotationYY< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_one_axis2.
  template <typename T>
  struct is_one_axis2_helper< qclab::qgates::RotationZZ< T > >
  : std::true_type { } ;

  /// Returns type of is_one_axis2.
  template <typename T>
  struct is_one_axis2
  : is_one_axis2_helper< T >::type { } ;

  /// Checks if T is a one axis rotation gate acting on two qubits.
  template <typename T>
  inline constexpr bool is_one_axis2_v = is_one_axis2< T >::value ;


  /// Checks if T is a one axis rotation gate.
  template <typename T>
  inline constexpr bool is_one_axis_v = ( is_one_axis1< T >::value ||
                                          is_one_axis2< T >::value ) ;


  /// Default helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper
  : std::false_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationX< T > ,
                              qclab::qgates::RotationX< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationX< T > ,
                              qclab::qgates::RotationXX< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationXX< T > ,
                              qclab::qgates::RotationX< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationXX< T > ,
                              qclab::qgates::RotationXX< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationY< T > ,
                              qclab::qgates::RotationY< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationY< T > ,
                              qclab::qgates::RotationYY< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationYY< T > ,
                              qclab::qgates::RotationY< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationYY< T > ,
                              qclab::qgates::RotationYY< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationZ< T > ,
                              qclab::qgates::RotationZ< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationZ< T > ,
                              qclab::qgates::RotationZZ< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationZZ< T > ,
                              qclab::qgates::RotationZ< U > >
  : std::true_type { } ;

  /// Template specialized helper class for is_same_axis.
  template <typename T, typename U>
  struct is_same_axis_helper< qclab::qgates::RotationZZ< T > ,
                              qclab::qgates::RotationZZ< U > >
  : std::true_type { } ;

  /// Returns type of is_same_axis.
  template <typename T, typename U>
  struct is_same_axis
  : is_same_axis_helper< T , U >::type { } ;

  /// Checks if T is a one axis rotation gate acting on two qubits.
  template <typename T, typename U>
  inline constexpr bool is_same_axis_v = is_same_axis< T , U >::value ;


  /// Default helper class for is_two_axes.
  template <typename T>
  struct is_two_axes_helper
  : std::false_type { } ;

  /// Template specialized helper class for is_two_axes.
  template <typename T>
  struct is_two_axes_helper< f3c::qgates::RotationXY< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_two_axes.
  template <typename T>
  struct is_two_axes_helper< f3c::qgates::RotationXZ< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_two_axes.
  template <typename T>
  struct is_two_axes_helper< f3c::qgates::RotationYZ< T > >
  : std::true_type { } ;

  /// Returns type of is_two_axes.
  template <typename T>
  struct is_two_axes
  : is_two_axes_helper< T >::type { } ;

  /// Checks if T is a two axes rotation gate.
  template <typename T>
  inline constexpr bool is_two_axes_v = is_two_axes< T >::value ;

// TODO: c++20
//  /// Concept: TwoAxesTurnoverable.
//  template <class T>
//  concept TwoAxesTurnoverable = is_two_axes_v< T > ;


  /// Default helper class for is_TF_two_axes.
  template <typename T>
  struct is_TF_two_axes_helper
  : std::false_type { } ;

  /// Template specialized helper class for is_TF_two_axes.
  template <typename T>
  struct is_TF_two_axes_helper< f3c::qgates::RotationTFXY< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_TF_two_axes.
  template <typename T>
  struct is_TF_two_axes_helper< f3c::qgates::RotationTFXZ< T > >
  : std::true_type { } ;

  /// Template specialized helper class for is_TF_two_axes.
  template <typename T>
  struct is_TF_two_axes_helper< f3c::qgates::RotationTFYZ< T > >
  : std::true_type { } ;

  /// Returns type of is_TF_two_axes.
  template <typename T>
  struct is_TF_two_axes
  : is_TF_two_axes_helper< T >::type { } ;

  /// Checks if T is a transverse fiel two axes rotation gate.
  template <typename T>
  inline constexpr bool is_TF_two_axes_v = is_TF_two_axes< T >::value ;

// TODO: c++20
//  /// Concept: TFTwoAxesTurnoverable.
//  template <class T>
//  concept TFTwoAxesTurnoverable = is_TF_two_axes_v< T > ;

} // namespace f3c

#endif

