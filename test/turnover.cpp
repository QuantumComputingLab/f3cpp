#include <gtest/gtest.h>
#include "f3c/turnover.hpp"
#include "qclab/QCircuit.hpp"
#include <random>

template <typename G1, typename G2>
void test_f3c_turnover_one_axis_1q() {

  using T = typename G1::value_type ;
  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // turnover VEE --> HAT
  //
  {
    G1 gate1( 0 , 0 ) ;
    G2 gate2( 0 , 0 ) ;
    G1 gate3( 0 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubit() , 0 ) ;
    EXPECT_EQ( gateB.qubit() , 0 ) ;
    EXPECT_EQ( gateC.qubit() , 0 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gate1 ) ) ;
    vee.push_back( std::make_unique< G2 >( gate2 ) ) ;
    vee.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gateA ) ) ;
    hat.push_back( std::make_unique< G1 >( gateB ) ) ;
    hat.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 10*eps ) ;
  }

  {
    G1 gate1( 1 , 0.1 ) ;
    G2 gate2( 1 , 0.2 ) ;
    G1 gate3( 1 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubit() , 1 ) ;
    EXPECT_EQ( gateB.qubit() , 1 ) ;
    EXPECT_EQ( gateC.qubit() , 1 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gate1 ) ) ;
    vee.push_back( std::make_unique< G2 >( gate2 ) ) ;
    vee.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gateA ) ) ;
    hat.push_back( std::make_unique< G1 >( gateB ) ) ;
    hat.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 10*eps ) ;
  }


  //
  // turnover HAT --> VEE
  //
  {
    G2 gate1( 1 , 0 ) ;
    G1 gate2( 1 , 0 ) ;
    G2 gate3( 1 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubit() , 1 ) ;
    EXPECT_EQ( gateB.qubit() , 1 ) ;
    EXPECT_EQ( gateC.qubit() , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gate1 ) ) ;
    hat.push_back( std::make_unique< G1 >( gate2 ) ) ;
    hat.push_back( std::make_unique< G2 >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gateA ) ) ;
    vee.push_back( std::make_unique< G2 >( gateB ) ) ;
    vee.push_back( std::make_unique< G1 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 10*eps ) ;
  }

  {
    G2 gate1( 0 , 0.1 ) ;
    G1 gate2( 0 , 0.2 ) ;
    G2 gate3( 0 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubit() , 0 ) ;
    EXPECT_EQ( gateB.qubit() , 0 ) ;
    EXPECT_EQ( gateC.qubit() , 0 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gate1 ) ) ;
    hat.push_back( std::make_unique< G1 >( gate2 ) ) ;
    hat.push_back( std::make_unique< G2 >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gateA ) ) ;
    vee.push_back( std::make_unique< G2 >( gateB ) ) ;
    vee.push_back( std::make_unique< G1 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 10*eps ) ;
  }

}


template <typename G1, typename G2>
void test_f3c_turnover_one_axis_2q() {

  using T = typename G1::value_type ;
  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  //           -- O -- X -- O --   -- X -- O -- X --
  // turnover     O         O    =         O
  //           -- O ------- O --   ------- O -------
  //
  if constexpr ( f3c::is_one_axis1_v< G2 > ) {
    G1 gate1( 0 , 1 , 0 ) ;
    G2 gate2( 0 , 0 ) ;
    G1 gate3( 0 , 1 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 1 ) ;
    EXPECT_EQ( gateB.nbQubits() , 2 ) ;
    EXPECT_EQ( gateC.nbQubits() , 1 ) ;
    EXPECT_EQ( gateA.qubit() , 0 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubit() , 0 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }

  if constexpr ( f3c::is_one_axis1_v< G2 > ) {
    G1 gate1( 0 , 1 , 0.1 ) ;
    G2 gate2( 0 , 0.2 ) ;
    G1 gate3( 0 , 1 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 1 ) ;
    EXPECT_EQ( gateB.nbQubits() , 2 ) ;
    EXPECT_EQ( gateC.nbQubits() , 1 ) ;
    EXPECT_EQ( gateA.qubit() , 0 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubit() , 0 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }


  //
  //           -- O ------- O --   ------- O -------
  // turnover     O         O    =         O
  //           -- O -- X -- O --   -- X -- O -- X --
  //
  if constexpr ( f3c::is_one_axis1_v< G2 > ) {
    G1 gate1( 0 , 1 , 0 ) ;
    G2 gate2( 1 , 0 ) ;
    G1 gate3( 0 , 1 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 1 ) ;
    EXPECT_EQ( gateB.nbQubits() , 2 ) ;
    EXPECT_EQ( gateC.nbQubits() , 1 ) ;
    EXPECT_EQ( gateA.qubit() , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubit() , 1 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }

  if constexpr ( f3c::is_one_axis1_v< G2 > ) {
    G1 gate1( 0 , 1 , 0.1 ) ;
    G2 gate2( 1 , 0.2 ) ;
    G1 gate3( 0 , 1 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 1 ) ;
    EXPECT_EQ( gateB.nbQubits() , 2 ) ;
    EXPECT_EQ( gateC.nbQubits() , 1 ) ;
    EXPECT_EQ( gateA.qubit() , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubit() , 1 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }


  //
  //           -- X -- O -- X --   -- O -- X -- O --
  // turnover          O         =    O         O
  //           ------- O -------   ---O-------- O --
  //
  if constexpr ( f3c::is_one_axis1_v< G1 > ) {
    G1 gate1( 0 , 0 ) ;
    G2 gate2( 0 , 1 , 0 ) ;
    G1 gate3( 0 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 2 ) ;
    EXPECT_EQ( gateB.nbQubits() , 1 ) ;
    EXPECT_EQ( gateC.nbQubits() , 2 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubit() , 0 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }

  if constexpr ( f3c::is_one_axis1_v< G1 > ) {
    G1 gate1( 0 , 0.1 ) ;
    G2 gate2( 0 , 1 , 0.2 ) ;
    G1 gate3( 0 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 2 ) ;
    EXPECT_EQ( gateB.nbQubits() , 1 ) ;
    EXPECT_EQ( gateC.nbQubits() , 2 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubit() , 0 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }


  //
  //           ------- O -------   -- O ------- O --
  // turnover          O         =    O         O
  //           -- X -- O -- X --   ---O--- X -- O --
  //
  if constexpr ( f3c::is_one_axis1_v< G1 > ) {
    G1 gate1( 1 , 0 ) ;
    G2 gate2( 0 , 1 , 0 ) ;
    G1 gate3( 1 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 2 ) ;
    EXPECT_EQ( gateB.nbQubits() , 1 ) ;
    EXPECT_EQ( gateC.nbQubits() , 2 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubit() , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }

  if constexpr ( f3c::is_one_axis1_v< G1 > ) {
    G1 gate1( 1 , 0.1 ) ;
    G2 gate2( 0 , 1 , 0.2 ) ;
    G1 gate3( 1 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.nbQubits() , 2 ) ;
    EXPECT_EQ( gateB.nbQubits() , 1 ) ;
    EXPECT_EQ( gateC.nbQubits() , 2 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubit() , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  left( 3 ) ;
    left.push_back( std::make_unique< G1 >( gate1 ) ) ;
    left.push_back( std::make_unique< G2 >( gate2 ) ) ;
    left.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  right( 3 ) ;
    right.push_back( std::make_unique< G2 >( gateA ) ) ;
    right.push_back( std::make_unique< G1 >( gateB ) ) ;
    right.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( left , right ) , 0.0 , 10*eps ) ;
  }

}


template <typename G1, typename G2>
void test_f3c_turnover_one_axis_3q() {

  using T = typename G1::value_type ;
  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // turnover VEE --> HAT
  //
  {
    G1 gate1( 0 , 1 , 0 ) ;
    G2 gate2( 1 , 2 , 0 ) ;
    G1 gate3( 0 , 1 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gate1 ) ) ;
    vee.push_back( std::make_unique< G2 >( gate2 ) ) ;
    vee.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gateA ) ) ;
    hat.push_back( std::make_unique< G1 >( gateB ) ) ;
    hat.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 10*eps ) ;
  }

  {
    G1 gate1( 0 , 1 , 0.1 ) ;
    G2 gate2( 1 , 2 , 0.2 ) ;
    G1 gate3( 0 , 1 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gate1 ) ) ;
    vee.push_back( std::make_unique< G2 >( gate2 ) ) ;
    vee.push_back( std::make_unique< G1 >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gateA ) ) ;
    hat.push_back( std::make_unique< G1 >( gateB ) ) ;
    hat.push_back( std::make_unique< G2 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 10*eps ) ;
  }


  //
  // turnover HAT --> VEE
  //
  {
    G2 gate1( 1 , 2 , 0 ) ;
    G1 gate2( 0 , 1 , 0 ) ;
    G2 gate3( 1 , 2 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gate1 ) ) ;
    hat.push_back( std::make_unique< G1 >( gate2 ) ) ;
    hat.push_back( std::make_unique< G2 >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gateA ) ) ;
    vee.push_back( std::make_unique< G2 >( gateB ) ) ;
    vee.push_back( std::make_unique< G1 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 10*eps ) ;
  }

  {
    G2 gate1( 1 , 2 , 0.1 ) ;
    G1 gate2( 0 , 1 , 0.2 ) ;
    G2 gate3( 1 , 2 , 0.3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G2 >( gate1 ) ) ;
    hat.push_back( std::make_unique< G1 >( gate2 ) ) ;
    hat.push_back( std::make_unique< G2 >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G1 >( gateA ) ) ;
    vee.push_back( std::make_unique< G2 >( gateB ) ) ;
    vee.push_back( std::make_unique< G1 >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 10*eps ) ;
  }

}


template <typename G>
void test_f3c_turnover_two_axes() {

  using T = typename G::value_type ;
  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // turnover VEE --> HAT
  //
  {
    G gate1( 0 , 1 , 0 , 0 ) ;
    G gate2( 1 , 2 , 0 , 0 ) ;
    G gate3( 0 , 1 , 0 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 10*eps ) ;
  }

  {
    G gate1( 0 , 1 , 0.1 , 0.6 ) ;
    G gate2( 1 , 2 , 0.2 , 0.7 ) ;
    G gate3( 0 , 1 , 0.3 , 0.8 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 10*eps ) ;
  }


  //
  // turnover HAT --> VEE
  //
  {
    G gate1( 1 , 2 , 0 , 0 ) ;
    G gate2( 0 , 1 , 0 , 0 ) ;
    G gate3( 1 , 2 , 0 , 0 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 10*eps ) ;
  }

  {
    G gate1( 1 , 2 , 0.1 , 0.6 ) ;
    G gate2( 0 , 1 , 0.2 , 0.7 ) ;
    G gate3( 1 , 2 , 0.3 , 0.8 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 10*eps ) ;
  }

}


template <typename G>
void test_f3c_turnover_TF() {

  using T = typename G::value_type ;
  using R = qclab::real_t< T > ;
  using TFXY = f3c::qgates::RotationTFXYMatrix< T > ;

  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // TFXY: turnover VEE --> HAT
  //
  {
    TFXY g1( 0 , 1 , 1 , 1 , 0 , 0 ) ;
    TFXY g2( 1 , 2 , 1 , 1 , 0 , 0 ) ;
    TFXY g3( 0 , 1 , 1 , 1 , 0 , 0 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 10*eps ) ;
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

    const T a3 = T(  7.573626678732417e-01 , -2.573742994005928e-01 ) ;
    const T b3 = T(  4.489593082886771e-01 ,  2.867411770124348e-01 ) ;
    const T c3 = T( -1.135746492799616e-01 , -8.386392764159148e-01 ) ;
    const T d3 = T( -5.811765837989824e-01 , -1.496463757119556e-01 ) ;

    TFXY g1( 0 , 1 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 1 , 2 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 0 , 1 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 100*eps ) ;
  }

  {
    const T a1 = 1 ;
    const T b1 = 1 ;
    const T c1 = 0 ;
    const T d1 = 0 ;

    const T a2 = T(  5.402590199290686e-01 , -8.414035709520180e-01 ) ;
    const T b2 = T(  9.999198857036123e-01 ,                      0 ) ;
    const T c2 = T(                      0 , -1.265788981128995e-02 ) ;
    const T d2 = T( -1.065124700509599e-02 , -6.839087052464791e-03 ) ;

    const T a3 = T( -1.265788981129948e-02 ,  3.105155021998485e-16 ) ;
    const T b3 = T( -1.265788981129948e-02 ,  3.105155021998485e-16 ) ;
    const T c3 = T( -8.414035709520178e-01 , -5.402590199290686e-01 ) ;
    const T d3 = T( -8.414035709520178e-01 , -5.402590199290686e-01 ) ;

    TFXY g1( 0 , 1 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 1 , 2 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 0 , 1 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 100*eps ) ;
  }

  {
    const T a1 = T(  9.800665778008439e-01 , -1.986693307868722e-01 ) ;
    const T b1 = T(  9.999999999587806e-01 ,  0.000000000000000e+00 ) ;
    const T c1 = T(  0.000000000000000e+00 , -9.079573740362129e-06 ) ;
    const T d1 = T( -1.803832838902155e-06 , -8.898586763973914e-06 ) ;

    const T a2 = T( -2.865803626278042e-05 ,  8.277844676001192e-06 ) ;
    const T b2 = T(  3.033605905429311e-05 ,  0.000000000000000e+00 ) ;
    const T c2 = T(  1.376119312981997e-06 ,  9.999999995389149e-01 ) ;
    const T d2 = T(  4.918273284784795e-06 , -9.999999995430024e-01 ) ;

    const T a3 = T(  9.210609937377550e-01 , -3.894183422524573e-01 ) ;
    const T b3 = T(  9.999999997443537e-01 ,  0.000000000000000e+00 ) ;
    const T c3 = T(  1.306322728480737e-05 ,  1.845656558177031e-05 ) ;
    const T d3 = T( -8.983471348101463e-06 ,  2.124769800254132e-05 ) ;

    TFXY g1( 0 , 1 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 1 , 2 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 0 , 1 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 100*eps ) ;
  }

  {
    const T a1 = T(  8.040304886458259e-08 , -5.749184594341116e-08 ) ;
    const T b1 = T(  7.373138834590875e-08 ,  6.583042601477045e-08 ) ;
    const T c1 = T(  9.983858291235649e-01 , -5.679556501384782e-02 ) ;
    const T d1 = T(  9.985685513389021e-01 , -5.348689817999129e-02 ) ;

    const T a2 = T( -5.169668546804159e-08 ,  8.424615913096427e-08 ) ;
    const T b2 = T(  9.595972106656034e-01 ,  2.813773147977348e-01 ) ;
    const T c2 = T(  4.051146504038741e-08 , -9.142658913183499e-08 ) ;
    const T d2 = T(  2.800326978703757e-01 ,  9.599904625169091e-01 ) ;

    const T a3 = T(  7.450242420396489e-01 , -6.670373893367871e-01 ) ;
    const T b3 = T(  1.602724700240165e-01 ,  9.870728115759195e-01 ) ;
    const T c3 = T(  5.914974679167009e-08 , -8.063068556375551e-08 ) ;
    const T d3 = T(  9.992877767346097e-08 ,  3.773511983273288e-09 ) ;

    TFXY g1( 0 , 1 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 1 , 2 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 0 , 1 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 100*eps ) ;
  }

  if constexpr ( std::is_same< qclab::real_t< T > , float >::value ) {
    // only float
    const T a1 = T(  8.27902615e-01 ,  5.60871899e-01 ) ;
    const T b1 = T(  0.00000000e+00 ,  0.00000000e+00 ) ;
    const T c1 = T(  9.35573339e-01 ,  3.53132486e-01 ) ;
    const T d1 = T(  9.74279146e-08 , -2.25344436e-08 ) ;

    const T a2 = T(  0.00000000e+00 , -0.00000000e+00 ) ;
    const T b2 = T(  4.95276242e-01 , -8.68735551e-01 ) ;
    const T c2 = T(  5.82292400e-08 ,  8.12979479e-08 ) ;
    const T d2 = T(  9.69858825e-01 , -2.43667617e-01 ) ;

    const T a3 = T(  7.63498485e-01 ,  6.45809590e-01 ) ;
    const T b3 = T(  8.77520084e-01 , -4.79539930e-01 ) ;
    const T c3 = T(  9.97455700e-08 ,  7.12884151e-09 ) ;
    const T d3 = T(  9.97505580e-08 , -7.05879932e-09 ) ;

    TFXY g1( 0 , 1 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 1 , 2 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 0 , 1 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 2 ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gate1 ) ) ;
    vee.push_back( std::make_unique< G >( gate2 ) ) ;
    vee.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gateA ) ) ;
    hat.push_back( std::make_unique< G >( gateB ) ) ;
    hat.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( vee , hat ) , 0.0 , 100*eps ) ;
  }


  //
  // TFXY: turnover HAT --> VEE
  //
  {
    TFXY g1( 1 , 2 , 1 , 1 , 0 , 0 ) ;
    TFXY g2( 0 , 1 , 1 , 1 , 0 , 0 ) ;
    TFXY g3( 1 , 2 , 1 , 1 , 0 , 0 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 10*eps ) ;
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

    const T a3 = T(  7.573626678732417e-01 , -2.573742994005928e-01 ) ;
    const T b3 = T(  4.489593082886771e-01 ,  2.867411770124348e-01 ) ;
    const T c3 = T( -1.135746492799616e-01 , -8.386392764159148e-01 ) ;
    const T d3 = T( -5.811765837989824e-01 , -1.496463757119556e-01 ) ;

    TFXY g1( 1 , 2 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 0 , 1 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 1 , 2 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 100*eps ) ;
  }

  {
    const T a1 = T(  8.508996043460090e-01 , -5.253283385881590e-01 ) ;
    const T b1 = T(  9.524673879925145e-01 , -3.046405665874241e-01 ) ;
    const T c1 = T(  9.822975476171381e-08 , -1.873273283462558e-08 ) ;
    const T d1 = T(  9.995811443211189e-08 , -2.894021281318401e-09 ) ;

    const T a2 = T(  9.750362944107126e-01 , -2.220455461877498e-01 ) ;
    const T b2 = T(  3.462360770782283e-01 , -9.381474185486358e-01 ) ;
    const T c2 = T(  1.887027753396164e-08 , -9.820342471518578e-08 ) ;
    const T d2 = T(  5.175058712180079e-08 ,  8.556796557444150e-08 ) ;

    const T a3 = T(  7.057856552777438e-08 ,  6.919991838682351e-08 ) ;
    const T b3 = T( -9.433603801253276e-09 ,  9.839191905854035e-08 ) ;
    const T c3 = T(  9.994792221685084e-01 , -3.226893945303590e-02 ) ;
    const T d3 = T(  3.044015169411465e-01 ,  9.525438134195818e-01 ) ;

    TFXY g1( 1 , 2 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 0 , 1 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 1 , 2 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 100*eps ) ;
  }

  {
    const T a1 = T(  5.597260981924069e-01 , -8.286776785954199e-01 ) ;
    const T b1 = T(  9.677577706900544e-01 ,  2.518827053789244e-01 ) ;
    const T c1 = T(  9.478166376834505e-08 ,  3.188159678096891e-08 ) ;
    const T d1 = T(  7.941370552610752e-08 , -6.077387082137010e-08 ) ;

    const T a2 = T(  7.278876416526158e-01 ,  6.856964205312605e-01 ) ;
    const T b2 = T(  6.746369175766853e-02 ,  9.977217298898607e-01 ) ;
    const T c2 = T(  2.761296116729983e-08 ,  9.611204074190232e-08 ) ;
    const T d2 = T(  9.661397083863748e-08 ,  2.580195029044911e-08 ) ;

    const T a3 = T( -5.034454048154690e-09 , -9.871482608574254e-08 ) ;
    const T b3 = T(  9.892381829124894e-01 ,  1.463141054990464e-01 ) ;
    const T c3 = T(  4.327581657942693e-08 ,  9.015100498266059e-08 ) ;
    const T d3 = T(  9.936461727682634e-01 , -1.125490263972667e-01 ) ;

    TFXY g1( 1 , 2 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 0 , 1 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 1 , 2 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 100*eps ) ;
  }

  {
    const T a1 = T( 6.340924154984683e-01 , -7.732572719395001e-01 ) ;
    const T b1 = T( 0.000000000000000e+00 ,  0.000000000000000e+00 ) ;
    const T c1 = T( 9.713077446841060e-01 , -2.378261237052726e-01 ) ;
    const T d1 = T( 9.999532192647693e-13 ,  9.672604717653607e-15 ) ;

    const T a2 = T( 0.000000000000000e+00 ,  0.000000000000000e+00 ) ;
    const T b2 = T( 0.000000000000000e+00 ,  0.000000000000000e+00 ) ;
    const T c2 = T( 9.471195219713482e-01 ,  3.208809921150908e-01 ) ;
    const T d2 = T( 8.690578371831881e-01 ,  4.947104967862310e-01 ) ;

    const T a3 = T( 0.000000000000000e+00 , -0.000000000000000e+00 ) ;
    const T b3 = T( 8.761394757430743e-01 , -4.820576926516690e-01 ) ;
    const T c3 = T( 9.788529481287986e-13 , -2.045651630643877e-13 ) ;
    const T d3 = T( 9.994141846336670e-01 , -3.422407855622338e-02 ) ;

    TFXY g1( 1 , 2 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 0 , 1 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 1 , 2 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 1000*eps ) ;
  }

  {
    const T a1 = T(  8.771666498303764e-01 , -4.801860768757814e-01 ) ;
    const T b1 = T(  0.000000000000000e+00 ,  0.000000000000000e+00 ) ;
    const T c1 = T(  1.120006420710291e-01 ,  9.937081343008506e-01 ) ;
    const T d1 = T(  5.033522803015200e-13 , -8.640812935802164e-13 ) ;

    const T a2 = T(  0.000000000000000e+00 ,  0.000000000000000e+00 ) ;
    const T b2 = T(  0.000000000000000e+00 , -0.000000000000000e+00 ) ;
    const T c2 = T(  7.827050908579476e-01 , -6.223927544123983e-01 ) ;
    const T d2 = T(  8.697745471024860e-01 ,  4.934493258812553e-01 ) ;

    const T a3 = T(  0.000000000000000e+00 ,  0.000000000000000e+00 ) ;
    const T b3 = T(  9.982129849640763e-01 ,  5.975647788406616e-02 ) ;
    const T c3 = T( -5.855328274366165e-14 ,  9.982842847004760e-13 ) ;
    const T d3 = T( -2.079633343912144e-01 , -9.781366221284673e-01 ) ;

    TFXY g1( 1 , 2 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 0 , 1 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 1 , 2 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 100*eps ) ;
  }

  if constexpr ( std::is_same< qclab::real_t< T > , float >::value ) {
    // only float
    const T a1 = T(  5.32970428e-01 , -8.46133887e-01 ) ;
    const T b1 = T(  5.81028819e-01 ,  8.13883006e-01 ) ;
    const T c1 = T(  9.99100939e-06 ,  4.23935091e-07 ) ;
    const T d1 = T(  9.04379339e-06 , -4.26729320e-06 ) ;

    const T a2 = T(  0.00000000e+00 , -0.00000000e+00 ) ;
    const T b2 = T(  7.86645710e-01 , -6.17404699e-01 ) ;
    const T c2 = T( -7.57578845e-06 , -6.52743619e-06 ) ;
    const T d2 = T(  9.24812197e-01 ,  3.80423963e-01 ) ;

    const T a3 = T(  5.49587726e-01 ,  8.35435986e-01 ) ;
    const T b3 = T(  9.84273254e-01 ,  1.76652744e-01 ) ;
    const T c3 = T(  5.84085501e-06 , -8.11692098e-06 ) ;
    const T d3 = T(  9.60437682e-06 , -2.78494894e-06 ) ;

    TFXY g1( 1 , 2 , a1 , b1 , c1 , d1 ) ;
    TFXY g2( 0 , 1 , a2 , b2 , c2 , d2 ) ;
    TFXY g3( 1 , 2 , a3 , b3 , c3 , d3 ) ;
    G gate1( g1 ) ;
    G gate2( g2 ) ;
    G gate3( g3 ) ;
    auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 , gate2 , gate3 ) ;
    EXPECT_EQ( gateA.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateA.qubits()[1] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[0] , 1 ) ;
    EXPECT_EQ( gateB.qubits()[1] , 2 ) ;
    EXPECT_EQ( gateC.qubits()[0] , 0 ) ;
    EXPECT_EQ( gateC.qubits()[1] , 1 ) ;

    qclab::QCircuit< T >  hat( 3 ) ;
    hat.push_back( std::make_unique< G >( gate1 ) ) ;
    hat.push_back( std::make_unique< G >( gate2 ) ) ;
    hat.push_back( std::make_unique< G >( gate3 ) ) ;

    qclab::QCircuit< T >  vee( 3 ) ;
    vee.push_back( std::make_unique< G >( gateA ) ) ;
    vee.push_back( std::make_unique< G >( gateB ) ) ;
    vee.push_back( std::make_unique< G >( gateC ) ) ;

    EXPECT_NEAR( qclab::nrmF( hat , vee ) , 0.0 , 100*eps ) ;
  }

}


template <typename R>
std::tuple< std::complex< R > ,
            std::complex< R > > test_f3c_turnover_genSU2( const R tol ) {

  using T = std::complex< R > ;

  std::random_device           rd ;
  std::mt19937                 gen( rd() ) ;
  std::normal_distribution<>   randn ;

  const R b = 1. / std::sqrt( 1. + tol*tol ) ;
  const R a = std::sqrt( 1. - b*b ) ;
  const R ra = randn( gen ) ;
  const R rb = randn( gen ) ;

  return { a * T( std::cos( ra ) , std::sin( ra ) ) ,
           b * T( std::cos( rb ) , std::sin( rb ) ) } ;

}


template <typename R>
void test_f3c_turnover_TFXY( const char type , const R min , const R max ) {

  using T = std::complex< R > ;
  using TFXY = f3c::qgates::RotationTFXYMatrix< T > ;

  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  const std::array< R , 2 >  ratios( { min , max } ) ;

  TFXY gate1 ;
  TFXY gate2 ;
  TFXY gate3 ;

  if ( type == 'v' ) {
    // vee --> hat
    gate1 = TFXY( 0 , 1 , 1 , 1 , 0 , 0 ) ;
    gate2 = TFXY( 1 , 2 , 1 , 1 , 0 , 0 ) ;
    gate3 = TFXY( 0 , 1 , 1 , 1 , 0 , 0 ) ;
  } else {
    // hat --> vee
    gate1 = TFXY( 1 , 2 , 1 , 1 , 0 , 0 ) ;
    gate2 = TFXY( 0 , 1 , 1 , 1 , 0 , 0 ) ;
    gate3 = TFXY( 1 , 2 , 1 , 1 , 0 , 0 ) ;
  }

  for ( int i1 = 0; i1 < 2; i1++ ) {
    for ( int i2 = 0; i2 < 2; i2++ ) {
      for ( int i3 = 0; i3 < 2; i3++ ) {
        for ( int i4 = 0; i4 < 2; i4++ ) {
          for ( int i5 = 0; i5 < 2; i5++ ) {
            for ( int i6 = 0; i6 < 2; i6++ ) {

              // update gates
              const auto [ a1 , d1 ] = test_f3c_turnover_genSU2( ratios[i1] ) ;
              const auto [ b1 , c1 ] = test_f3c_turnover_genSU2( ratios[i2] ) ;
              gate1.update( a1 , b1 , c1 , d1 ) ;
              const auto [ a2 , d2 ] = test_f3c_turnover_genSU2( ratios[i3] ) ;
              const auto [ b2 , c2 ] = test_f3c_turnover_genSU2( ratios[i4] ) ;
              gate2.update( a2 , b2 , c2 , d2 ) ;
              const auto [ a3 , d3 ] = test_f3c_turnover_genSU2( ratios[i5] ) ;
              const auto [ b3 , c3 ] = test_f3c_turnover_genSU2( ratios[i6] ) ;
              gate3.update( a3 , b3 , c3 , d3 ) ;

              // turnover
              const auto [ gateA , gateB , gateC ] = f3c::turnover( gate1 ,
                                                                    gate2 ,
                                                                    gate3 ) ;

              // quantum circuit in
              qclab::QCircuit< T >  cin( 3 ) ;
              cin.push_back( std::make_unique< TFXY >( gate1 ) ) ;
              cin.push_back( std::make_unique< TFXY >( gate2 ) ) ;
              cin.push_back( std::make_unique< TFXY >( gate3 ) ) ;

              // quantum circuit out
              qclab::QCircuit< T >  cout( 3 ) ;
              cout.push_back( std::make_unique< TFXY >( gateA ) ) ;
              cout.push_back( std::make_unique< TFXY >( gateB ) ) ;
              cout.push_back( std::make_unique< TFXY >( gateC ) ) ;

              // Frobenius norm
              const auto nrmF = qclab::nrmF( cin , cout ) ;
              if ( nrmF > 100*eps ) {
                std::printf( "type = %c: (%e,%e)\n" , type , min , max ) ;
                std::printf( "  a1 = %22.15e + 1i * %22.15e\n" , std::real( a1 ) , std::imag( a1 ) ) ;
                std::printf( "  b1 = %22.15e + 1i * %22.15e\n" , std::real( b1 ) , std::imag( b1 ) ) ;
                std::printf( "  c1 = %22.15e + 1i * %22.15e\n" , std::real( c1 ) , std::imag( c1 ) ) ;
                std::printf( "  d1 = %22.15e + 1i * %22.15e\n" , std::real( d1 ) , std::imag( d1 ) ) ;
                std::printf( "  a2 = %22.15e + 1i * %22.15e\n" , std::real( a2 ) , std::imag( a2 ) ) ;
                std::printf( "  b2 = %22.15e + 1i * %22.15e\n" , std::real( b2 ) , std::imag( b2 ) ) ;
                std::printf( "  c2 = %22.15e + 1i * %22.15e\n" , std::real( c2 ) , std::imag( c2 ) ) ;
                std::printf( "  d2 = %22.15e + 1i * %22.15e\n" , std::real( d2 ) , std::imag( d2 ) ) ;
                std::printf( "  a3 = %22.15e + 1i * %22.15e\n" , std::real( a3 ) , std::imag( a3 ) ) ;
                std::printf( "  b3 = %22.15e + 1i * %22.15e\n" , std::real( b3 ) , std::imag( b3 ) ) ;
                std::printf( "  c3 = %22.15e + 1i * %22.15e\n" , std::real( c3 ) , std::imag( c3 ) ) ;
                std::printf( "  d3 = %22.15e + 1i * %22.15e\n" , std::real( d3 ) , std::imag( d3 ) ) ;
                std::cout << std::endl ;
                std::printf( "  A1 = %22.15e + 1i * %22.15e\n" , std::real( gateA.a() ) , std::imag( gateA.a() ) ) ;
                std::printf( "  B1 = %22.15e + 1i * %22.15e\n" , std::real( gateA.b() ) , std::imag( gateA.b() ) ) ;
                std::printf( "  C1 = %22.15e + 1i * %22.15e\n" , std::real( gateA.c() ) , std::imag( gateA.c() ) ) ;
                std::printf( "  D1 = %22.15e + 1i * %22.15e\n" , std::real( gateA.d() ) , std::imag( gateA.d() ) ) ;
                std::printf( "  A2 = %22.15e + 1i * %22.15e\n" , std::real( gateB.a() ) , std::imag( gateB.a() ) ) ;
                std::printf( "  B2 = %22.15e + 1i * %22.15e\n" , std::real( gateB.b() ) , std::imag( gateB.b() ) ) ;
                std::printf( "  C2 = %22.15e + 1i * %22.15e\n" , std::real( gateB.c() ) , std::imag( gateB.c() ) ) ;
                std::printf( "  D2 = %22.15e + 1i * %22.15e\n" , std::real( gateB.d() ) , std::imag( gateB.d() ) ) ;
                std::printf( "  A3 = %22.15e + 1i * %22.15e\n" , std::real( gateC.a() ) , std::imag( gateC.a() ) ) ;
                std::printf( "  B3 = %22.15e + 1i * %22.15e\n" , std::real( gateC.b() ) , std::imag( gateC.b() ) ) ;
                std::printf( "  C3 = %22.15e + 1i * %22.15e\n" , std::real( gateC.c() ) , std::imag( gateC.c() ) ) ;
                std::printf( "  D3 = %22.15e + 1i * %22.15e\n" , std::real( gateC.d() ) , std::imag( gateC.d() ) ) ;
              }
              EXPECT_NEAR( qclab::nrmF( cin , cout ) , 0.0 , 100*eps ) ;

            }
          }
        }
      }
    }
  }

}


template <typename T>
void test_f3c_turnover() {

  test_f3c_turnover_one_axis_1q< qclab::qgates::RotationX< T > ,
                                 qclab::qgates::RotationY< T > >() ;
  test_f3c_turnover_one_axis_1q< qclab::qgates::RotationY< T > ,
                                 qclab::qgates::RotationX< T > >() ;
  test_f3c_turnover_one_axis_1q< qclab::qgates::RotationX< T > ,
                                 qclab::qgates::RotationZ< T > >() ;
  test_f3c_turnover_one_axis_1q< qclab::qgates::RotationZ< T > ,
                                 qclab::qgates::RotationX< T > >() ;
  test_f3c_turnover_one_axis_1q< qclab::qgates::RotationY< T > ,
                                 qclab::qgates::RotationZ< T > >() ;
  test_f3c_turnover_one_axis_1q< qclab::qgates::RotationZ< T > ,
                                 qclab::qgates::RotationY< T > >() ;

  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationX< T > ,
                                 qclab::qgates::RotationYY< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationX< T > ,
                                 qclab::qgates::RotationZZ< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationY< T > ,
                                 qclab::qgates::RotationXX< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationY< T > ,
                                 qclab::qgates::RotationZZ< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationZ< T > ,
                                 qclab::qgates::RotationXX< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationZ< T > ,
                                 qclab::qgates::RotationYY< T > >() ;

  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationXX< T > ,
                                 qclab::qgates::RotationY< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationXX< T > ,
                                 qclab::qgates::RotationZ< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationYY< T > ,
                                 qclab::qgates::RotationX< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationYY< T > ,
                                 qclab::qgates::RotationZ< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationZZ< T > ,
                                 qclab::qgates::RotationX< T > >() ;
  test_f3c_turnover_one_axis_2q< qclab::qgates::RotationZZ< T > ,
                                 qclab::qgates::RotationY< T > >() ;

  test_f3c_turnover_one_axis_3q< qclab::qgates::RotationXX< T > ,
                                 qclab::qgates::RotationYY< T > >() ;
  test_f3c_turnover_one_axis_3q< qclab::qgates::RotationYY< T > ,
                                 qclab::qgates::RotationXX< T > >() ;
  test_f3c_turnover_one_axis_3q< qclab::qgates::RotationXX< T > ,
                                 qclab::qgates::RotationZZ< T > >() ;
  test_f3c_turnover_one_axis_3q< qclab::qgates::RotationZZ< T > ,
                                 qclab::qgates::RotationXX< T > >() ;
  test_f3c_turnover_one_axis_3q< qclab::qgates::RotationYY< T > ,
                                 qclab::qgates::RotationZZ< T > >() ;
  test_f3c_turnover_one_axis_3q< qclab::qgates::RotationZZ< T > ,
                                 qclab::qgates::RotationYY< T > >() ;

  test_f3c_turnover_two_axes< f3c::qgates::RotationXY< T > >() ;
  test_f3c_turnover_two_axes< f3c::qgates::RotationXZ< T > >() ;
  test_f3c_turnover_two_axes< f3c::qgates::RotationYZ< T > >() ;

  test_f3c_turnover_TF< f3c::qgates::RotationTFXY< T > >() ;
  test_f3c_turnover_TF< f3c::qgates::RotationTFXZ< T > >() ;
  test_f3c_turnover_TF< f3c::qgates::RotationTFYZ< T > >() ;
  test_f3c_turnover_TF< f3c::qgates::RotationTFXYMatrix< T > >() ;

}


/*
 * float
 */
TEST( f3c_turnover , complex_float ) {
  test_f3c_turnover< std::complex< float > >() ;
  test_f3c_turnover_TFXY< float >( 'v' , 1e-3 , 1e3 ) ;
  test_f3c_turnover_TFXY< float >( 'h' , 1e-3 , 1e3 ) ;
  test_f3c_turnover_TFXY< float >( 'v' , 1e-5 , 1e5 ) ;
  test_f3c_turnover_TFXY< float >( 'h' , 1e-5 , 1e5 ) ;
}

/*
 * double
 */
TEST( f3c_turnover , complex_double ) {
  test_f3c_turnover< std::complex< double > >() ;
  test_f3c_turnover_TFXY< double >( 'v' , 1e-3 , 1e3 ) ;
  test_f3c_turnover_TFXY< double >( 'h' , 1e-3 , 1e3 ) ;
  test_f3c_turnover_TFXY< double >( 'v' , 1e-7 , 1e7 ) ;
  test_f3c_turnover_TFXY< double >( 'h' , 1e-7 , 1e7 ) ;
  test_f3c_turnover_TFXY< double >( 'v' , 1e-12 , 1e12 ) ;
  test_f3c_turnover_TFXY< double >( 'h' , 1e-12 , 1e12 ) ;
}

