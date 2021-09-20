//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_RotationTFXZ_hpp
#define f3c_qgates_RotationTFXZ_hpp

#include "f3c/qgates/TFTwoAxesQRotationGate2.hpp"
#include "f3c/qasm.hpp"
#include "qclab/qgates/RotationX.hpp"
#include "qclab/qgates/RotationY.hpp"
#include "qclab/qgates/RotationZ.hpp"
#include "qclab/qgates/CNOT.hpp"
#include "qclab/QCircuit.hpp"

namespace f3c {

  namespace qgates {

    /**
     * \class RotationTFXZ
     * \brief 2-qubit transverse field rotation gate about XY.
     */
    template <typename T>
    class RotationTFXZ : public TFTwoAxesQRotationGate2< T >
    {

      public:
        /// Real value type of this TFXZ-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this TFXZ-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a TFXZ-rotation gate on
         *        qubits 0 and 1 with parameters \f$\theta_0 = \theta_1 =
         *        \theta_2 = \theta_3 = \theta_4 = \theta_5 = 0\f$.
         */
        RotationTFXZ()
        : TFTwoAxesQRotationGate2< T >()
        { } // RotationTFXZ()

        /**
         * \brief Constructs a TFXZ-rotation gate on qubits 0 and 1 with
         *        the given quantum rotations
         *          `rot0` = \f$\theta_0\f$, `rot1` = \f$\theta_1\f$,
         *          `rot2` = \f$\theta_2\f$, `rot3` = \f$\theta_3\f$,
         *          `rot4` = \f$\theta_4\f$, `rot5` = \f$\theta_5\f$.
         */
        RotationTFXZ( const rotation_type& rot0 , const rotation_type& rot1 ,
                      const rotation_type& rot2 , const rotation_type& rot3 ,
                      const rotation_type& rot4 , const rotation_type& rot5 )
        : TFTwoAxesQRotationGate2< T >( 0 , 1 , rot0 , rot1 , rot2 , rot3 ,
                                        rot4 , rot5)
        { } // RotationTFXZ(rot0,rot1,rot2,rot3,rot4,rot5)

        /**
         * \brief Constructs a TFXZ-rotation gate on qubits 0 and 1 with
         *        the given values
         *          `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *          `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *          `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        RotationTFXZ( const real_type theta0 , const real_type theta1 ,
                      const real_type theta2 , const real_type theta3 ,
                      const real_type theta4 , const real_type theta5 )
        : TFTwoAxesQRotationGate2< T >( 0 , 1 ,
                                        rotation_type( theta0 ) ,
                                        rotation_type( theta1 ) ,
                                        rotation_type( theta2 ) ,
                                        rotation_type( theta3 ) ,
                                        rotation_type( theta4 ) ,
                                        rotation_type( theta5 ) )
        { } // RotationTFXZ(theta0,theta1,theta2,theta3,theta4,theta5)

        /**
         * \brief Constructs a TFXZ-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with quantum rotations
         *          `rot0` = \f$\theta_0\f$, `rot1` = \f$\theta_1\f$,
         *          `rot2` = \f$\theta_2\f$, `rot3` = \f$\theta_3\f$,
         *          `rot4` = \f$\theta_4\f$, `rot5` = \f$\theta_5\f$.
         */
        RotationTFXZ( const int qubit0 , const int qubit1 ,
                      const rotation_type& rot0 , const rotation_type& rot1 ,
                      const rotation_type& rot2 , const rotation_type& rot3 ,
                      const rotation_type& rot4 , const rotation_type& rot5 )
        : TFTwoAxesQRotationGate2< T >( qubit0 , qubit1 , rot0 , rot1 ,
                                        rot2 , rot3 , rot4 , rot5 )
        { } // RotationTFXZ(qubit0,qubit1,rot0,rot1,rot2,rot3,rot4,rot5)

        /**
         * \brief Constructs a TFXZ-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with values
         *          `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *          `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *          `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        RotationTFXZ( const int qubit0 , const int qubit1 ,
                      const real_type theta0 , const real_type theta1 ,
                      const real_type theta2 , const real_type theta3 ,
                      const real_type theta4 , const real_type theta5 )
        : TFTwoAxesQRotationGate2< T >( qubit0 , qubit1 ,
                                        rotation_type( theta0 ) ,
                                        rotation_type( theta1 ) ,
                                        rotation_type( theta2 ) ,
                                        rotation_type( theta3 ) ,
                                        rotation_type( theta4 ) ,
                                        rotation_type( theta5 ) )
        { } // RotationTFXZ(qubit0,qubit1,theta0,theta1,theta2,theta3,theta4,theta5)

        /**
         * \brief Constructs a TFXZ-rotation gate from the given TFXY-rotation
         *        matrix gate `gate`.
         */
        RotationTFXZ( const RotationTFXYMatrix< T >& gate )
        : TFTwoAxesQRotationGate2< T >( gate )
        { } // RotationTFXZ(gate)

        // nbQubits

        // fixed

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        /// Returns the unitary matrix corresponding to this TFXZ-rotation gate.
        qclab::dense::SquareMatrix< T > matrix() const override {
          using RX = qclab::qgates::RotationX< T > ;
          using RY = qclab::qgates::RotationY< T > ;
          using RZ = qclab::qgates::RotationZ< T > ;
          using CNOT = qclab::qgates::CNOT< T > ;
          qclab::QCircuit< T > circuit( 2 ) ;
          circuit.push_back( std::make_unique< RY >( 0 , this->theta0() ) ) ;
          circuit.push_back( std::make_unique< RY >( 1 , this->theta1() ) ) ;
          circuit.push_back( std::make_unique< CNOT >( 0 , 1 ) ) ;
          circuit.push_back( std::make_unique< RX >( 0 , this->theta2() ) ) ;
          circuit.push_back( std::make_unique< RZ >( 1 , this->theta3() ) ) ;
          circuit.push_back( std::make_unique< CNOT >( 0 , 1 ) ) ;
          circuit.push_back( std::make_unique< RY >( 0 , this->theta4() ) ) ;
          circuit.push_back( std::make_unique< RY >( 1 , this->theta5() ) ) ;
          return circuit.matrix() ;
        }

        // apply

        // print

        /**
         * \brief Writes the QASM code of this TFXZ-rotation gate to the given
         *        `stream`.
         */
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmTFRxz( this->qubits_[0] + offset ,
                               this->qubits_[1] + offset ,
                               this->theta0() , this->theta1() ,
                               this->theta2() , this->theta3() ,
                               this->theta4() , this->theta5() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        /// Checks if `other` equals this TFXZ-rotation gate.
        inline bool equals( const qclab::QObject< T >& other ) const override {
          using TFXZ = RotationTFXZ< T > ;
          if ( const TFXZ* p = dynamic_cast< const TFXZ* >( &other ) ) {
            const auto [ other0 , other1 , other2 ,
                         other3 , other4 , other5 ] = p->rotations() ;
            const auto [ this0 , this1 , this2 ,
                         this3 , this4 , this5 ] = this->rotations() ;
            return ( other0 == this0 ) && ( other1 == this1 ) &&
                   ( other2 == this2 ) && ( other3 == this3 ) &&
                   ( other4 == this4 ) && ( other5 == this5 ) ;
          }
          return false ;
        }

        // rotations

        // thetas

        // update(rot0,rot1,rot2,rot3,rot4,rot5)

        // update(theta0,theta1,theta2,theta3,theta4,theta5)

        // a

        // b

        // c

        // d

        /// Multiplies `rhs` to this 2-qubit TFXZ-rotation gate.
        inline RotationTFXZ< T >& operator*=( const RotationTFXZ< T >& rhs ) {
          assert( this->qubits()[0] == rhs.qubits()[0] ) ;
          assert( this->qubits()[1] == rhs.qubits()[1] ) ;
          using TFXYMatrix = RotationTFXYMatrix< T > ;
          // fuse
          RotationTFXYMatrix< T > matrix( *this ) ;
          matrix *= RotationTFXYMatrix< T >( rhs ) ;
          const RotationTFXZ< T > tmp( matrix ) ;
          // update
          this->rot_[0] = tmp.rotation0() ;
          this->rot_[1] = tmp.rotation1() ;
          this->rot_[2] = tmp.rotation2() ;
          this->rot_[3] = tmp.rotation3() ;
          this->rot_[4] = tmp.rotation4() ;
          this->rot_[5] = tmp.rotation5() ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationTFXZ< T > operator*( RotationTFXZ< T > lhs ,
                                            const RotationTFXZ< T >& rhs ) {
          assert( lhs.qubits()[0] == rhs.qubits()[0] ) ;
          assert( lhs.qubits()[1] == rhs.qubits()[1] ) ;
          lhs *= rhs ;
          return lhs ;
        }

    } ; // class RotationTFXZ

  } // namespace qgates

} // namespace f3c

#endif

