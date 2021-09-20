//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_RotationYZ_hpp
#define f3c_qgates_RotationYZ_hpp

#include "f3c/qgates/TwoAxesQRotationGate2.hpp"
#include "f3c/qasm.hpp"

namespace f3c {

  namespace qgates {

    /**
     * \class RotationYZ
     * \brief 2-qubit rotation gate about YZ.
     */
    template <typename T>
    class RotationYZ : public TwoAxesQRotationGate2< T >
    {

      public:
        /// Real value type of this YZ-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this YZ-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a YZ-rotation gate on
         *        qubits 0 and 1 with parameters \f$\theta_0 = \theta_1 = 0\f$.
         */
        RotationYZ()
        : TwoAxesQRotationGate2< T >()
        { } // RotationYZ()

        /**
         * \brief Constructs a YZ-rotation gate on qubits 0 and 1 with
         *        the given quantum rotations `rot0` = \f$\theta_0\f$ and
         *                                    `rot1` = \f$\theta_1\f$.
         */
        RotationYZ( const rotation_type& rot0 , const rotation_type& rot1 )
        : TwoAxesQRotationGate2< T >( 0 , 1 , rot0 , rot1 )
        { } // RotationYZ(rot0,rot1)

        /**
         * \brief Constructs a YZ-rotation gate on qubits 0 and 1 with
         *        the given values `theta0` = \f$\theta_0\f$ and
         *                         `theta1` = \f$\theta_1\f$.
         */
        RotationYZ( const real_type theta0 , const real_type theta1 )
        : TwoAxesQRotationGate2< T >( 0 , 1 , theta0 , theta1 )
        { } // RotationYZ(theta0,theta1)

        /**
         * \brief Constructs a YZ-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with quantum rotations `rot0` = \f$\theta_0\f$
         *        and `rot1` = \f$\theta_1\f$.
         */
        RotationYZ( const int qubit0 , const int qubit1 ,
                    const rotation_type& rot0 , const rotation_type& rot1 )
        : TwoAxesQRotationGate2< T >( qubit0 , qubit1 , rot0 , rot1 )
        { } // RotationYZ(qubit0,qubit1,rot0,rot1)

        /**
         * \brief Constructs a YZ-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with values `theta0` = \f$\theta_0\f$
         *        and `theta1` = \f$\theta_1\f$.
         */
        RotationYZ( const int qubit0 , const int qubit1 ,
                    const real_type theta0 , const real_type theta1 )
        : TwoAxesQRotationGate2< T >( qubit0 , qubit1 , theta0 , theta1 )
        { } // RotationYZ(qubit0,qubit1,theta0,theta1)

        // nbQubits

        // fixed

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        /// Returns the unitary matrix corresponding to this YZ-rotation gate.
        qclab::dense::SquareMatrix< T > matrix() const override {
          const auto [ rot0 , rot1 ] = this->rotations() ;
          const T cp = T(  rot1.cos() ,  rot1.sin() ) * rot0.cos() ;
          const T cm = T(  rot1.cos() , -rot1.sin() ) * rot0.cos() ;
          const T sp = T(  rot1.sin() , -rot1.cos() ) * rot0.sin() ;
          const T sm = T(  rot1.sin() ,  rot1.cos() ) * rot0.sin() ;
          return qclab::dense::SquareMatrix< T >( cm ,  0 ,  0 , sm ,
                                                   0 , cp , sp ,  0 ,
                                                   0 , sp , cp ,  0 ,
                                                  sm ,  0 ,  0 , cm ) ;
        }

        // apply

        // print

        /// Writes the QASM code of this YZ-rotation gate to the given `stream`.
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          auto qubits = this->qubits() ;
          auto [ theta0 , theta1 ] = this->thetas() ;
          stream << qasmRyz( qubits[0] + offset , qubits[1] + offset ,
                             theta0 , theta1 ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        /// Checks if `other` equals this YZ-rotation gate.
        inline bool equals( const qclab::QObject< T >& other ) const override {
          using YZ = RotationYZ< T > ;
          if ( const YZ* p = dynamic_cast< const YZ* >( &other ) ) {
            const auto [ other0 , other1 ] = p->rotations() ;
            const auto [  this0 ,  this1 ] = this->rotations() ;
            return ( other0 == this0 ) && ( other1 == this1 ) ;
          }
          return false ;
        }

        // rotations

        // thetas

        // update(rot0,rot1)

        // update(theta0,theta1)

        /// Multiplies `rhs` to this 2-qubit YZ-rotation gate.
        inline RotationYZ< T >& operator*=( const RotationYZ< T >& rhs ) {
          assert( this->qubits()[0] == rhs.qubits()[0] ) ;
          assert( this->qubits()[1] == rhs.qubits()[1] ) ;
          this->rotations_[0] *= rhs.rotations_[0] ;
          this->rotations_[1] *= rhs.rotations_[1] ;
          return *this ;
        }

        /// Multiplies the inverse of `rhs` to this 2-qubit YZ-rotation gate.
        inline RotationYZ< T >& operator/=( const RotationYZ< T >& rhs ) {
          assert( this->qubits()[0] == rhs.qubits()[0] ) ;
          assert( this->qubits()[1] == rhs.qubits()[1] ) ;
          this->rotations_[0] /= rhs.rotations_[0] ;
          this->rotations_[1] /= rhs.rotations_[1] ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationYZ< T > operator*( RotationYZ< T > lhs ,
                                          const RotationYZ< T >& rhs ) {
          assert( lhs.qubits()[0] == rhs.qubits()[0] ) ;
          assert( lhs.qubits()[1] == rhs.qubits()[1] ) ;
          lhs *= rhs ;
          return lhs ;
        }

        /// Multiplies `lhs` and the inverse of `rhs`.
        friend RotationYZ< T > operator/( RotationYZ< T > lhs ,
                                          const RotationYZ< T >& rhs ) {
          assert( lhs.qubits()[0] == rhs.qubits()[0] ) ;
          assert( lhs.qubits()[1] == rhs.qubits()[1] ) ;
          lhs /= rhs ;
          return lhs ;
        }

        /// Returns the inverse this 2-qubit YZ-rotation gate.
        inline RotationYZ< T > inv() const {
          RotationYZ< T > rotation( -this->rotations_[0].angle() ,
                                    -this->rotations_[1].angle() ) ;
          return rotation ;
        }

    } ; // class RotationYZ

  } // namespace qgates

} // namespace f3c

#endif

