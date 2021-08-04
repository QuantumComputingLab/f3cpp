//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_RotationTFXY_hpp
#define f3c_qgates_RotationTFXY_hpp

#include "qclab/qgates/QGate2.hpp"
#include "qclab/QRotation.hpp"
#include "f3c/turnoverSU2.hpp"
#include "f3c/qasm.hpp"
#include <array>

namespace f3c {

  namespace qgates {

    template <typename T>
    class RotationTFXYMatrix ;

    /**
     * \class RotationTFXY
     * \brief 2-qubit transverse field rotation gate about XY.
     */
    template <typename T>
    class RotationTFXY : public qclab::qgates::QGate2< T >
    {

      public:
        /// Real value type of this TFXY-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this TFXY-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a TFXY-rotation gate on
         *        qubits 0 and 1 with parameters \f$\theta_0 = \theta_1 =
         *        \theta_2 = \theta_3 = \theta_4 = \theta_5 = 0\f$.
         */
        RotationTFXY()
        : qubits_( { 0 , 1 } )
        { } // RotationTFXY()

        /**
         * \brief Constructs a TFXY-rotation gate on qubits 0 and 1 with
         *        the given quantum rotations
         *   `rot0` = \f$\theta_0+\theta_1\f$, `rot1` = \f$\theta_0-\theta_1\f$,
         *   `rot2` = \f$\theta_2+\theta_3\f$, `rot3` = \f$\theta_2-\theta_3\f$,
         *   `rot4` = \f$\theta_4+\theta_5\f$, `rot5` = \f$\theta_4-\theta_5\f$.
         */
        RotationTFXY( const rotation_type& rot0 , const rotation_type& rot1 ,
                      const rotation_type& rot2 , const rotation_type& rot3 ,
                      const rotation_type& rot4 , const rotation_type& rot5 )
        : qubits_( { 0 , 1 } )
        , rot_( { rot0 , rot1 , rot2 , rot3 , rot4 , rot5 } )
        { } // RotationTFXY(rot0,rot1,rot2,rot3,rot4,rot5)

        /**
         * \brief Constructs a TFXY-rotation gate on qubits 0 and 1 with
         *        the given values
         *          `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *          `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *          `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        RotationTFXY( const real_type theta0 , const real_type theta1 ,
                      const real_type theta2 , const real_type theta3 ,
                      const real_type theta4 , const real_type theta5 )
        : qubits_( { 0 , 1 } )
        , rot_( { rotation_type( theta0 + theta1 ) ,
                  rotation_type( theta0 - theta1 ) ,
                  rotation_type( theta2 + theta3 ) ,
                  rotation_type( theta2 - theta3 ) ,
                  rotation_type( theta4 + theta5 ) ,
                  rotation_type( theta4 - theta5 ) } )
        { } // RotationTFXY(theta0,theta1,theta2,theta3,theta4,theta5)

        /**
         * \brief Constructs a TFXY-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with quantum rotations
         *   `rot0` = \f$\theta_0+\theta_1\f$, `rot1` = \f$\theta_0-\theta_1\f$,
         *   `rot2` = \f$\theta_2+\theta_3\f$, `rot3` = \f$\theta_2-\theta_3\f$,
         *   `rot4` = \f$\theta_4+\theta_5\f$, `rot5` = \f$\theta_4-\theta_5\f$.
         */
        RotationTFXY( const int qubit0 , const int qubit1 ,
                      const rotation_type& rot0 , const rotation_type& rot1 ,
                      const rotation_type& rot2 , const rotation_type& rot3 ,
                      const rotation_type& rot4 , const rotation_type& rot5 )
        : rot_( { rot0 , rot1 , rot2 , rot3 , rot4 , rot5 } )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // RotationTFXY(qubit0,qubit1,rot0,rot1,rot2,rot3,rot4,rot5)

        /**
         * \brief Constructs a TFXY-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with values
         *          `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *          `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *          `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        RotationTFXY( const int qubit0 , const int qubit1 ,
                      const real_type theta0 , const real_type theta1 ,
                      const real_type theta2 , const real_type theta3 ,
                      const real_type theta4 , const real_type theta5 )
        : rot_( { rotation_type( theta0 + theta1 ) ,
                  rotation_type( theta0 - theta1 ) ,
                  rotation_type( theta2 + theta3 ) ,
                  rotation_type( theta2 - theta3 ) ,
                  rotation_type( theta4 + theta5 ) ,
                  rotation_type( theta4 - theta5 ) } )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // RotationTFXY(qubit0,qubit1,theta0,theta1,theta2,theta3,theta4,theta5)

        /**
         * \brief Constructs a TFXY-rotation gate from the given TFXY-rotation
         *        gate `gate`.
         */
        RotationTFXY( const RotationTFXYMatrix< T >& gate )
        {
          // qubits
          setQubits( &(gate.qubits()[0]) ) ;
          // rotations
          const T a = gate.a() ;
          const T b = gate.b() ;
          const T c = gate.c() ;
          const T d = gate.d() ;
          using R = real_type ;
          const R anglea = std::atan2( std::imag(a) , std::real(a) ) ;
          const R angleb = std::atan2( std::imag(b) , std::real(b) ) ;
          const R anglec = std::atan2( std::imag(c) , std::real(c) ) ;
          const R angled = std::atan2( std::imag(d) , std::real(d) ) ;
          const R pi2 = 2 * std::atan(1) ;
          rot_[0] = rotation_type( -anglea - angled - pi2 ) ;
          rot_[1] = rotation_type( -angleb - anglec - pi2 ) ;
          rot_[2] = rotation_type( 2*std::atan2( std::abs(c) , std::abs(b) ) ) ;
          rot_[3] = rotation_type( 2*std::atan2( std::abs(d) , std::abs(a) ) ) ;
          rot_[4] = rotation_type( -anglea + angled + pi2 ) ;
          rot_[5] = rotation_type( -angleb + anglec + pi2 ) ;
        } // RotationTFXY(gate)

        // nbQubits

        /// Checks if this TFXY-rotation gate is fixed.
        inline bool fixed() const override { return false ; }

        /// Checks if this TFXY-rotation gate is controlled.
        inline bool controlled() const override { return false ; }

        /// Returns the first qubit of this TFXY-rotation gate.
        inline int qubit() const override { return qubits_[0] ; }

        // setQubit

        /// Returns the qubits of this TFXY-rotation gate in ascending order.
        std::vector< int > qubits() const override {
          return std::vector< int >( { qubits_[0] , qubits_[1] } ) ;
        }

        /// Sets the qubits of this TFXY-rotation gate.
        inline void setQubits( const int* qubits ) override {
          assert( qubits[0] >= 0 ) ; assert( qubits[1] >= 0 ) ;
          assert( qubits[0] != qubits[1] ) ;
          qubits_[0] = std::min( qubits[0] , qubits[1] ) ;
          qubits_[1] = std::max( qubits[0] , qubits[1] ) ;
        }

        /// Returns the unitary matrix corresponding to this TFXY-rotation gate.
        qclab::dense::SquareMatrix< T > matrix() const override {
          const T a = this->a() ;
          const T b = this->b() ;
          const T c = this->c() ;
          const T d = this->d() ;
          using M = qclab::dense::SquareMatrix< T > ;
          return M( a , 0 ,       0       , -std::conj(d) ,
                    0 , b , -std::conj(c) ,       0       ,
                    0 , c ,  std::conj(b) ,       0       ,
                    d , 0 ,       0       ,  std::conj(a) ) ;
        }

        // apply

        // print

        /**
         * \brief Writes the QASM code of this TFXY-rotation gate to the given
         *        `stream`.
         */
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmTFRxy( qubits_[0] + offset , qubits_[1] + offset ,
                               theta0() , theta1() , theta2() ,
                               theta3() , theta4() , theta5() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        /// Checks if `other` equals this TFXY-rotation gate.
        inline bool equals( const qclab::QObject< T >& other ) const override {
          using TFXY = RotationTFXY< T > ;
          if ( const TFXY* p = dynamic_cast< const TFXY* >( &other ) ) {
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

        /**
         * \brief Returns the quantum rotation \f$\theta_0 + \theta_1\f$ of this
         *        TFXY-rotation gate.
         */
        inline const rotation_type& rotation0() const { return rot_[0] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_0 - \theta_1\f$ of this
         *        TFXY-rotation gate.
         */
        inline const rotation_type& rotation1() const { return rot_[1] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_2 + \theta_3\f$ of this
         *        TFXY-rotation gate.
         */
        inline const rotation_type& rotation2() const { return rot_[2] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_2 - \theta_3\f$ of this
         *        TFXY-rotation gate.
         */
        inline const rotation_type& rotation3() const { return rot_[3] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_4 + \theta_5\f$ of this
         *        TFXY-rotation gate.
         */
        inline const rotation_type& rotation4() const { return rot_[4] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_4 - \theta_5\f$ of this
         *        TFXY-rotation gate.
         */
        inline const rotation_type& rotation5() const { return rot_[5] ; } ;

        /**
         * \brief Returns the quantum rotations
         *          \f$\theta_0 + \theta_1\f$, \f$\theta_0 - \theta_1\f$,
         *          \f$\theta_2 + \theta_3\f$, \f$\theta_2 - \theta_3\f$,
         *          \f$\theta_4 + \theta_5\f$, \f$\theta_4 - \theta_5\f$,
         *        of this TFXY-rotation gate.
         */
        std::tuple< const rotation_type& ,
                    const rotation_type& ,
                    const rotation_type& ,
                    const rotation_type& ,
                    const rotation_type& ,
                    const rotation_type& > rotations() const {
          return { rot_[0] , rot_[1] , rot_[2] , rot_[3] , rot_[4] , rot_[5] } ;
        }

        /**
         * \brief Returns the numerical value \f$\theta_0\f$ of this
         *        TFXY-rotation gate.
         */
        inline real_type theta0() const {
          return ( rot_[0].theta() + rot_[1].theta() ) / 2 ;
        }

        /**
         * \brief Returns the numerical value \f$\theta_1\f$ of this
         *        TFXY-rotation gate.
         */
        inline real_type theta1() const {
          return ( rot_[0].theta() - rot_[1].theta() ) / 2 ;
        }

        /**
         * \brief Returns the numerical value \f$\theta_2\f$ of this
         *        TFXY-rotation gate.
         */
        inline real_type theta2() const {
          return ( rot_[2].theta() + rot_[3].theta() ) / 2 ;
        }

        /**
         * \brief Returns the numerical value \f$\theta_3\f$ of this
         *        TFXY-rotation gate.
         */
        inline real_type theta3() const {
          return ( rot_[2].theta() - rot_[3].theta() ) / 2 ;
        }

        /**
         * \brief Returns the numerical value \f$\theta_4\f$ of this
         *        TFXY-rotation gate.
         */
        inline real_type theta4() const {
          return ( rot_[4].theta() + rot_[5].theta() ) / 2 ;
        }

        /**
         * \brief Returns the numerical value \f$\theta_5\f$ of this
         *        TFXY-rotation gate.
         */
        inline real_type theta5() const {
          return ( rot_[4].theta() - rot_[5].theta() ) / 2 ;
        }

        /**
         * \brief Returns the numerical values \f$\theta_0\f$, \f$\theta_1\f$,
         *        \f$\theta_2\f$, \f$\theta_3\f$, \f$\theta_4\f$, and
         *        \f$\theta_5\f$ of this TFXY-rotation gate.
         */
        std::tuple< real_type , real_type , real_type ,
                    real_type , real_type , real_type > thetas() const {
          return { theta0() , theta1() , theta2() ,
                   theta3() , theta4() , theta5() } ;
        }

        /**
         * \brief Updates this TFXY-rotation gate with the given quantum
         *        rotations
         *   `rot0` = \f$\theta_0+\theta_1\f$, `rot1` = \f$\theta_0-\theta_1\f$,
         *   `rot2` = \f$\theta_2+\theta_3\f$, `rot3` = \f$\theta_2-\theta_3\f$,
         *   `rot4` = \f$\theta_4+\theta_5\f$, `rot5` = \f$\theta_4-\theta_5\f$.
         */
        void update( const rotation_type& rot0 , const rotation_type& rot1 ,
                     const rotation_type& rot2 , const rotation_type& rot3 ,
                     const rotation_type& rot4 , const rotation_type& rot5 ) {
          rot_[0] = rot0 ;
          rot_[1] = rot1 ;
          rot_[2] = rot2 ;
          rot_[3] = rot3 ;
          rot_[4] = rot4 ;
          rot_[5] = rot5 ;
        }

        /**
         * \brief Updates this TFXY-rotation gate with the given values
         *          `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *          `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *          `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        void update( const real_type theta0 , const real_type theta1 ,
                     const real_type theta2 , const real_type theta3 ,
                     const real_type theta4 , const real_type theta5 ) {
          rot_[0].update( theta0 + theta1 ) ;
          rot_[1].update( theta0 - theta1 ) ;
          rot_[2].update( theta2 + theta3 ) ;
          rot_[3].update( theta2 - theta3 ) ;
          rot_[4].update( theta4 + theta5 ) ;
          rot_[5].update( theta4 - theta5 ) ;
        }

        /// Returns the numerical value \f$a\f$ of this TFXY-rotation gate.
        inline T a() const {
          const auto tmp = rot_[0] * rot_[4] ;
          return T( rot_[3].cos() * tmp.cos() , -rot_[3].cos() * tmp.sin() ) ;
        }

        /// Returns the numerical value \f$b\f$ of this TFXY-rotation gate.
        inline T b() const {
          const auto tmp = rot_[1] * rot_[5] ;
          return T( rot_[2].cos() * tmp.cos() , -rot_[2].cos() * tmp.sin() ) ;
        }

        /// Returns the numerical value \f$c\f$ of this TFXY-rotation gate.
        inline T c() const {
          const auto tmp = rot_[1] / rot_[5] ;
          return T( -rot_[2].sin() * tmp.sin() , -rot_[2].sin() * tmp.cos() ) ;
        }

        /// Returns the numerical value \f$d\f$ of this TFXY-rotation gate.
        inline T d() const {
          const auto tmp = rot_[0] / rot_[4] ;
          return T( -rot_[3].sin() * tmp.sin() , -rot_[3].sin() * tmp.cos() ) ;
        }


        /// Multiplies `rhs` to this 2-qubit TFXY-rotation gate.
        inline RotationTFXY< T >& operator*=( const RotationTFXY< T >& rhs ) {
          assert( this->qubits()[0] == rhs.qubits()[0] ) ;
          assert( this->qubits()[1] == rhs.qubits()[1] ) ;
          // fuse
          rot_[4] *= rhs.rot_[0] ;
          rot_[5] *= rhs.rot_[1] ;
          // refactor the Euler decompositions for two SU(2)s
          const auto [ left1 , mid1 , right1 ]
                              = turnoverSU2( rot_[3] , rot_[4] , rhs.rot_[3] ) ;
          const auto [ left2 , mid2 , right2 ]
                              = turnoverSU2( rot_[2] , rot_[5] , rhs.rot_[2] ) ;
          // update
          rot_[0] *= left1 ;
          rot_[1] *= left2 ;
          rot_[2] = mid2 ;
          rot_[3] = mid1 ;
          rot_[4] = rhs.rot_[4] * right1 ;
          rot_[5] = rhs.rot_[5] * right2 ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationTFXY< T > operator*( RotationTFXY< T > lhs ,
                                            const RotationTFXY< T >& rhs ) {
          assert( lhs.qubits()[0] == rhs.qubits()[0] ) ;
          assert( lhs.qubits()[1] == rhs.qubits()[1] ) ;
          lhs *= rhs ;
          return lhs ;
        }

      protected:
        /// Qubits of this TFXY-rotation gate.
        std::array< int , 2 >            qubits_ ;
        /// Rotations of this TFXY-rotation gate.
        std::array< rotation_type , 6 >  rot_ ;

    } ; // class RotationTFXY

  } // namespace qgates

} // namespace f3c

#endif

