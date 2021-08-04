//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_RotationXY_hpp
#define f3c_qgates_RotationXY_hpp

#include "f3c/qgates/TwoAxesQRotationGate2.hpp"
#include "f3c/qasm.hpp"

namespace f3c {

  namespace qgates {

    /**
     * \class RotationXY
     * \brief 2-qubit rotation gate about XY.
     */
    template <typename T>
    class RotationXY : public TwoAxesQRotationGate2< T >
    {

      public:
        /// Real value type of this XY-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this XY-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs an XY-rotation gate on
         *        qubits 0 and 1 with parameters \f$\theta_0 = \theta_1 = 0\f$.
         */
        RotationXY()
        : TwoAxesQRotationGate2< T >()
        { } // RotationXY()

        /**
         * \brief Constructs an XY-rotation gate on qubits 0 and 1 with
         *        the given quantum rotations `rot0` = \f$\theta_0\f$ and
         *                                    `rot1` = \f$\theta_1\f$.
         */
        RotationXY( const rotation_type& rot0 , const rotation_type& rot1 )
        : TwoAxesQRotationGate2< T >( 0 , 1 , rot0 , rot1 )
        { } // RotationXY(rot0,rot1)

        /**
         * \brief Constructs an XY-rotation gate on qubits 0 and 1 with
         *        the given values `theta0` = \f$\theta_0\f$ and
         *                         `theta1` = \f$\theta_1\f$.
         */
        RotationXY( const real_type theta0 , const real_type theta1 )
        : TwoAxesQRotationGate2< T >( 0 , 1 , theta0 , theta1 )
        { } // RotationXY(theta0,theta1)

        /**
         * \brief Constructs an XY-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with quantum rotations `rot0` = \f$\theta_0\f$
         *        and `rot1` = \f$\theta_1\f$.
         */
        RotationXY( const int qubit0 , const int qubit1 ,
                    const rotation_type& rot0 , const rotation_type& rot1 )
        : TwoAxesQRotationGate2< T >( qubit0 , qubit1 , rot0 , rot1 )
        { } // RotationXY(qubit0,qubit1,rot0,rot1)

        /**
         * \brief Constructs an XY-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with values `theta0` = \f$\theta_0\f$
         *        and `theta1` = \f$\theta_1\f$.
         */
        RotationXY( const int qubit0 , const int qubit1 ,
                    const real_type theta0 , const real_type theta1 )
        : TwoAxesQRotationGate2< T >( qubit0 , qubit1 , theta0 , theta1 )
        { } // RotationXY(qubit0,qubit1,theta0,theta1)

        // nbQubits

        // fixed

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        /// Returns the unitary matrix corresponding to this XY-rotation gate.
        qclab::dense::SquareMatrix< T > matrix() const override {
          const auto [ rot0 , rot1 ] = this->rotations() ;
          const auto angle0 = rot0.angle() ;
          const auto angle1 = rot1.angle() ;
          const auto plus  = angle0 + angle1 ;
          const auto minus = angle0 - angle1 ;
          const T cp = plus.cos() ;
          const T sp = T(0,-plus.sin()) ;
          const T cm = minus.cos() ;
          const T sm = T(0,-minus.sin()) ;
          return qclab::dense::SquareMatrix< T >( cm ,  0 ,  0 , sm ,
                                                   0 , cp , sp ,  0 ,
                                                   0 , sp , cp ,  0 ,
                                                  sm ,  0 ,  0 , cm ) ;
        }

        // apply

        // print

        /// Writes the QASM code of this XY-rotation gate to the given `stream`.
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          auto qubits = this->qubits() ;
          auto [ theta0 , theta1 ] = this->thetas() ;
          stream << qasmRxy( qubits[0] + offset , qubits[1] + offset ,
                             theta0 , theta1 ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        /// Checks if `other` equals this XY-rotation gate.
        inline bool equals( const qclab::QObject< T >& other ) const override {
          using XY = RotationXY< T > ;
          if ( const XY* p = dynamic_cast< const XY* >( &other ) ) {
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

        /// Multiplies `rhs` to this 2-qubit XY-rotation gate.
        inline RotationXY< T >& operator*=( const RotationXY< T >& rhs ) {
          assert( this->qubits()[0] == rhs.qubits()[0] ) ;
          assert( this->qubits()[1] == rhs.qubits()[1] ) ;
          this->rotations_[0] *= rhs.rotations_[0] ;
          this->rotations_[1] *= rhs.rotations_[1] ;
          return *this ;
        }

        /// Multiplies the inverse of `rhs` to this 2-qubit XY-rotation gate.
        inline RotationXY< T >& operator/=( const RotationXY< T >& rhs ) {
          assert( this->qubits()[0] == rhs.qubits()[0] ) ;
          assert( this->qubits()[1] == rhs.qubits()[1] ) ;
          this->rotations_[0] /= rhs.rotations_[0] ;
          this->rotations_[1] /= rhs.rotations_[1] ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationXY< T > operator*( RotationXY< T > lhs ,
                                          const RotationXY< T >& rhs ) {
          assert( lhs.qubits()[0] == rhs.qubits()[0] ) ;
          assert( lhs.qubits()[1] == rhs.qubits()[1] ) ;
          lhs *= rhs ;
          return lhs ;
        }

        /// Multiplies `lhs` and the inverse of `rhs`.
        friend RotationXY< T > operator/( RotationXY< T > lhs ,
                                          const RotationXY< T >& rhs ) {
          assert( lhs.qubits()[0] == rhs.qubits()[0] ) ;
          assert( lhs.qubits()[1] == rhs.qubits()[1] ) ;
          lhs /= rhs ;
          return lhs ;
        }

        /// Returns the inverse this 2-qubit XY-rotation gate.
        inline RotationXY< T > inv() const {
          RotationXY< T > rotation( -this->rotations_[0].angle() ,
                                    -this->rotations_[1].angle() ) ;
          return rotation ;
        }

    } ; // class RotationXY

  } // namespace qgates

} // namespace f3c

#endif

