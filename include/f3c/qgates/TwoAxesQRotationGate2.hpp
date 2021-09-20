//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_TwoAxesQRotationGate2_hpp
#define f3c_qgates_TwoAxesQRotationGate2_hpp

#include "qclab/qgates/QGate2.hpp"
#include "qclab/QRotation.hpp"
#include <array>

namespace f3c {

  namespace qgates {

    /**
     * \class TwoAxesQRotationGate2
     * \brief Base class for 2-qubit rotation gates with 2-axes interactions.
     */
    template <typename T>
    class TwoAxesQRotationGate2 : public qclab::qgates::QGate2< T >
    {

      public:
        /// Real value type of this 2-axes rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this 2-axes rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a 2-axes rotation gate on
         *        qubits 0 and 1 with parameters \f$\theta_0 = \theta_1 = 0\f$.
         */
        TwoAxesQRotationGate2()
        : qubits_( { 0 , 1 } )
        , rotations_( { rotation_type() , rotation_type() } )
        { } // TwoAxesQRotationGate2()

        /**
         * \brief Constructs a 2-axes rotation gate on the given qubits `qubit0`
         *        and `qubit1` with quantum rotations `rot0` = \f$\theta_0\f$
         *        and `rot1` = \f$\theta_1\f$.
         */
        TwoAxesQRotationGate2( const int qubit0 , const int qubit1 ,
                               const rotation_type& rot0 ,
                               const rotation_type& rot1 )
        : rotations_( { rot0 , rot1 } )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // TwoAxesQRotationGate2(qubit0,qubit1,rot0,rot1)

        /**
         * \brief Constructs a 2-axes rotation gate on the given qubits `qubit0`
         *        and `qubit1` with quantum thetas `theta0` = \f$\theta_0\f$
         *        and `theta1` = \f$\theta_1\f$.
         */
        TwoAxesQRotationGate2( const int qubit0 , const int qubit1 ,
                               const real_type theta0 , const real_type theta1 )
        : rotations_( { rotation_type( theta0 ) , rotation_type( theta1 ) } )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // TwoAxesQRotationGate2(qubit0,qubit1,theta0,theta1)

        // nbQubits

        /// Checks if this 2-axes rotation gate is fixed.
        inline bool fixed() const override { return false ; }

        /// Checks if this 2-axes rotation gate is controlled.
        inline bool controlled() const override { return false ; }

        /// Returns the first qubit of this 2-axes rotation gate.
        inline int qubit() const override { return qubits_[0] ; }

        // setQubit

        /// Returns the qubits of this 2-axes rotation gate in ascending order.
        std::vector< int > qubits() const override {
          return std::vector< int >( { qubits_[0] , qubits_[1] } ) ;
        }

        /// Sets the qubits of this 2-axes rotation gate.
        inline void setQubits( const int* qubits ) override {
          assert( qubits[0] >= 0 ) ; assert( qubits[1] >= 0 ) ;
          assert( qubits[0] != qubits[1] ) ;
          qubits_[0] = std::min( qubits[0] , qubits[1] ) ;
          qubits_[1] = std::max( qubits[0] , qubits[1] ) ;
        }

        // matrix

        // apply

        // print

        // toQASM

        // operator==

        // operator!=

        // equals

        /**
         * \brief Returns the quantum rotations \f$\theta_0\f$ and
         *        \f$\theta_1\f$ of this 2-axes rotation gate.
         */
        inline std::tuple< const rotation_type& ,
                           const rotation_type& > rotations() const {
          return { rotations_[0] , rotations_[1] } ;
        }

        /**
         * \brief Returns the numerical values \f$\theta_0\f$ and \f$\theta_1\f$
         *        of this 2-axes rotation gate.
         */
        inline std::tuple< real_type , real_type > thetas() const {
          return { rotations_[0].theta() , rotations_[1].theta() } ;
        }

        /**
         * \brief Updates this 2-axes rotation gate with the given quantum
         *        rotations `rot0` = \f$\theta_0\f$ and `rot1` = \f$\theta_1\f$.
         */
        void update( const rotation_type& rot0 , const rotation_type& rot1 ) {
          rotations_[0] = rot0 ;
          rotations_[1] = rot1 ;
        }

        /**
         * \brief Updates this 2-axes rotation gate with the given values
         *        `theta0` = \f$\theta_0\f$ and `theta1` = \f$\theta_1\f$.
         */
        void update( const real_type theta0 , const real_type theta1 ) {
          rotations_[0].update( theta0 ) ;
          rotations_[1].update( theta1 ) ;
        }

      protected:
        /// Qubits of this 2-axes rotation gate.
        std::array< int , 2 >            qubits_ ;
        /// Quantum rotations of this 2-axes rotation gate.
        std::array< rotation_type , 2 >  rotations_ ;

    } ; // class TwoAxesQRotationGate2

  } // namespace qgates

} // namespace f3c

#endif

