//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_TFTwoAxesQRotationGate2_hpp
#define f3c_qgates_TFTwoAxesQRotationGate2_hpp

#include "qclab/qgates/QGate2.hpp"
#include "qclab/QRotation.hpp"
#include "f3c/qasm.hpp"
#include <array>

namespace f3c {

  namespace qgates {

    template <typename T>
    class RotationTFXYMatrix ;

    /**
     * \class TFTwoAxesQRotationGate2
     * \brief Base class for 2-qubit transverse field rotation gates with
     *        2-axes interactions.
     */
    template <typename T>
    class TFTwoAxesQRotationGate2 : public qclab::qgates::QGate2< T >
    {

      public:
        /// Real value type of this transverse field 2-axes rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this transverse field 2-axes rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a transverse field 2-axes
         *        rotation gate on qubits 0 and 1 with parameters \f$\theta_0 =
         *        \theta_1 = \theta_2 = \theta_3 = \theta_4 = \theta_5 = 0\f$.
         */
        TFTwoAxesQRotationGate2()
        : qubits_( { 0 , 1 } )
        { } // TFTwoAxesQRotationGate2()

        /**
         * \brief Constructs a transverse field 2-axes rotation gate on the
         *        given qubits `qubit0` and `qubit1` with quantum rotations
         *          `rot0` = \f$\theta_0\f$, `rot1` = \f$\theta_1\f$,
         *          `rot2` = \f$\theta_2\f$, `rot3` = \f$\theta_3\f$,
         *          `rot4` = \f$\theta_4\f$, `rot5` = \f$\theta_5\f$.
         */
        TFTwoAxesQRotationGate2( const int qubit0 , const int qubit1 ,
                                 const rotation_type& rot0 ,
                                 const rotation_type& rot1 ,
                                 const rotation_type& rot2 ,
                                 const rotation_type& rot3 ,
                                 const rotation_type& rot4 ,
                                 const rotation_type& rot5 )
        : rot_( { rot0 , rot1 , rot2 , rot3 , rot4 , rot5 } )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // TFTwoAxesQRotationGate2(qubit0,qubit1,rot0,rot1,rot2,rot3,rot4,rot5)

        /**
         * \brief Constructs a transverse field 2-axes rotation gate on the
         *        given qubits `qubit0` and `qubit1` with values
         *          `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *          `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *          `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        TFTwoAxesQRotationGate2( const int qubit0 , const int qubit1 ,
                                 const real_type theta0 ,
                                 const real_type theta1 ,
                                 const real_type theta2 ,
                                 const real_type theta3 ,
                                 const real_type theta4 ,
                                 const real_type theta5 )
        : rot_( { rotation_type( theta0 ) , rotation_type( theta1 ) ,
                  rotation_type( theta2 ) , rotation_type( theta3 ) ,
                  rotation_type( theta4 ) , rotation_type( theta5 ) } )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // TFTwoAxesQRotationGate2(qubit0,qubit1,theta0,theta1,theta2,theta3,theta4,theta5)

        /**
         * \brief Constructs a transverse field 2-axes rotation gate from the
         *        given transverse field XY-rotation matrix gate `gate`.
         */
        TFTwoAxesQRotationGate2( const RotationTFXYMatrix< T >& gate )
        {
          // qubits
          setQubits( &(gate.qubits()[0]) ) ;
          // rotations
          const T a = gate.a() ;
          const T b = gate.b() ;
          const T c = gate.c() ;
          const T d = gate.d() ;
          using R = real_type ;
          const R arga = std::arg( a ) ;
          const R argb = std::arg( b ) ;
          const R argc = std::arg( c ) ;
          const R argd = std::arg( d ) ;
          const R pi = 4 * std::atan(1) ;
          const R theta1 = std::atan2( std::abs(c) , std::abs(b) ) ;
          const R theta2 = std::atan2( std::abs(d) , std::abs(a) ) ;
          rot_[0] = rotation_type( ( -arga - argb - argc - argd - pi ) / 2 ) ;
          rot_[1] = rotation_type( ( -arga + argb + argc - argd      ) / 2 ) ;
          rot_[2] = rotation_type( theta1 + theta2 ) ;
          rot_[3] = rotation_type( theta1 - theta2 ) ;
          rot_[4] = rotation_type( ( -arga - argb + argc + argd + pi ) / 2 ) ;
          rot_[5] = rotation_type( ( -arga + argb - argc + argd      ) / 2 ) ;
        } // TFTwoAxesQRotationGate2(gate)

        // nbQubits

        /// Checks if this transverse field 2-axes rotation gate is fixed.
        inline bool fixed() const override { return false ; }

        /// Checks if this transverse field 2-axes rotation gate is controlled.
        inline bool controlled() const override { return false ; }

        /**
         * \brief Returns the first qubit of this transverse field 2-axes
         *        rotation gate.
         */
        inline int qubit() const override { return qubits_[0] ; }

        // setQubit

        /**
         * \brief Returns the qubits of this transverse field 2-axes rotation
         *        gate in ascending order.
         */
        std::vector< int > qubits() const override {
          return std::vector< int >( { qubits_[0] , qubits_[1] } ) ;
        }

        /// Sets the qubits of this transverse field 2-axes rotation gate.
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
         * \brief Returns the quantum rotation \f$\theta_0\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline const rotation_type& rotation0() const { return rot_[0] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_1\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline const rotation_type& rotation1() const { return rot_[1] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_2\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline const rotation_type& rotation2() const { return rot_[2] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_3\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline const rotation_type& rotation3() const { return rot_[3] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_4\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline const rotation_type& rotation4() const { return rot_[4] ; } ;

        /**
         * \brief Returns the quantum rotation \f$\theta_5\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline const rotation_type& rotation5() const { return rot_[5] ; } ;

        /**
         * \brief Returns the quantum rotations \f$\theta_0\f$, \f$\theta_1\f$,
         *        \f$\theta_2\f$, \f$\theta_3\f$, \f$\theta_4\f$, and
         *        \f$\theta_5\f$ of this transverse field 2-axes rotation gate.
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
         *        transverse field 2-axes rotation gate.
         */
        inline real_type theta0() const { return rot_[0].theta() ; }

        /**
         * \brief Returns the numerical value \f$\theta_1\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline real_type theta1() const { return rot_[1].theta() ; }

        /**
         * \brief Returns the numerical value \f$\theta_2\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline real_type theta2() const { return rot_[2].theta() ; }

        /**
         * \brief Returns the numerical value \f$\theta_3\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline real_type theta3() const { return rot_[3].theta() ; }

        /**
         * \brief Returns the numerical value \f$\theta_4\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline real_type theta4() const { return rot_[4].theta() ; }

        /**
         * \brief Returns the numerical value \f$\theta_5\f$ of this
         *        transverse field 2-axes rotation gate.
         */
        inline real_type theta5() const { return rot_[5].theta() ; }

        /**
         * \brief Returns the numerical values \f$\theta_0\f$, \f$\theta_1\f$,
         *        \f$\theta_2\f$, \f$\theta_3\f$, \f$\theta_4\f$, and
         *        \f$\theta_5\f$ of this transverse field 2-axes rotation gate.
         */
        std::tuple< real_type , real_type , real_type ,
                    real_type , real_type , real_type > thetas() const {
          return { theta0() , theta1() , theta2() ,
                   theta3() , theta4() , theta5() } ;
        }

        /**
         * \brief Updates this transverse field rotation gate with the given
         *        quantum rotations
         *          `rot0` = \f$\theta_0\f$, `rot1` = \f$\theta_1\f$,
         *          `rot2` = \f$\theta_2\f$, `rot3` = \f$\theta_3\f$,
         *          `rot4` = \f$\theta_4\f$, `rot5` = \f$\theta_5\f$.
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
         * \brief Updates this transverse field rotation gate with the given
         *        values `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *               `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *               `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        void update( const real_type theta0 , const real_type theta1 ,
                     const real_type theta2 , const real_type theta3 ,
                     const real_type theta4 , const real_type theta5 ) {
          rot_[0].update( theta0 ) ;
          rot_[1].update( theta1 ) ;
          rot_[2].update( theta2 ) ;
          rot_[3].update( theta3 ) ;
          rot_[4].update( theta4 ) ;
          rot_[5].update( theta5 ) ;
        }

        /**
         * \brief Returns the numerical value \f$a\f$ of this transverse field
         *        rotation gate.
         */
        inline T a() const {
          const auto scal = ( rot_[2] / rot_[3] ).cos() ;
          const auto theta = ( rot_[0] * rot_[1] ) * ( rot_[4] * rot_[5] ) ;
          return T( scal * theta.cos() , -scal * theta.sin() ) ;
        }

        /**
         * \brief Returns the numerical value \f$b\f$ of this transverse field
         *        rotation gate.
         */
        inline T b() const {
          const auto scal = ( rot_[2] * rot_[3] ).cos() ;
          const auto theta = ( rot_[0] / rot_[1] ) * ( rot_[4] / rot_[5] ) ;
          return T( scal * theta.cos() , -scal * theta.sin() ) ;
        }

        /**
         * \brief Returns the numerical value \f$c\f$ of this transverse field
         *        rotation gate.
         */
        inline T c() const {
          const auto scal = ( rot_[2] * rot_[3] ).sin() ;
          const auto theta = ( rot_[0] / rot_[1] ) / ( rot_[4] / rot_[5] ) ;
          return T( -scal * theta.sin() , -scal * theta.cos() ) ;
        }

        /**
         * \brief Returns the numerical value \f$d\f$ of this transverse field
         *        rotation gate.
         */
        inline T d() const {
          const auto scal = ( rot_[2] / rot_[3] ).sin() ;
          const auto theta = ( rot_[0] * rot_[1] ) / ( rot_[4] * rot_[5] ) ;
          return T( -scal * theta.sin() , -scal * theta.cos() ) ;
        }

      protected:
        /// Qubits of this transverse field rotation gate.
        std::array< int , 2 >            qubits_ ;
        /// Rotations of this transverse field rotation gate.
        std::array< rotation_type , 6 >  rot_ ;

    } ; // class TFTwoAxesQRotationGate2

  } // namespace qgates

} // namespace f3c

#endif

