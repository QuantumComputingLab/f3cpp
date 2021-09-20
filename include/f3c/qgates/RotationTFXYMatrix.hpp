//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qgates_RotationTFXYMatrix_hpp
#define f3c_qgates_RotationTFXYMatrix_hpp

#include "qclab/qgates/QGate2.hpp"
#include <array>

namespace f3c {

  namespace qgates {

    template <typename T>
    class RotationTFXY ;

    template <typename T>
    class RotationTFXZ ;

    template <typename T>
    class RotationTFYZ ;

    /**
     * \class RotationTFXYMatrix
     * \brief 2-qubit transverse field rotation gate about XY.
     */
    template <typename T>
    class RotationTFXYMatrix : public qclab::qgates::QGate2< T >
    {

      public:
        /// Real value type of this TFXY-rotation gate.
        using real_type = qclab::real_t< T > ;

        /**
         * \brief Default constructor. Constructs a TFXY-rotation gate on
         *        qubits 0 and 1 with parameters \f$a = b = 1\f$ and
         *        \f$c = d = 0\f$.
         */
        RotationTFXYMatrix()
        : qubits_( { 0 , 1 } )
        , v_( { 1 , 1 , 0 , 0 } )
        { } // RotationTFXYMatrix()

        /**
         * \brief Constructs a TFXY-rotation gate on qubits 0 and 1 with
         *        the given numerical values `a`, `b`, `c`, and `d`.
         */
        RotationTFXYMatrix( const T a , const T b , const T c , const T d )
        : qubits_( { 0 , 1 } )
        {
          update( a , b , c , d ) ;
        } // RotationTFXYMatrix(a,b,c,d)

        /**
         * \brief Constructs a TFXY-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with numerical values `a`, `b`, `c`, and `d`.
         */
        RotationTFXYMatrix( const int qubit0 , const int qubit1 ,
                            const T a , const T b , const T c , const T d )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
          update( a , b , c , d ) ;
        } // RotationTFXYMatrix(qubit0,qubit1,a,b,c,d)

        /**
         * \brief Constructs a TFXY-rotation gate on the given qubits `qubit0`
         *        and `qubit1` with values
         *          `theta0` = \f$\theta_0\f$, `theta1` = \f$\theta_1\f$,
         *          `theta2` = \f$\theta_2\f$, `theta3` = \f$\theta_3\f$,
         *          `theta4` = \f$\theta_4\f$, `theta5` = \f$\theta_5\f$.
         */
        RotationTFXYMatrix( const int qubit0 , const int qubit1 ,
                            const real_type theta0 , const real_type theta1 ,
                            const real_type theta2 , const real_type theta3 ,
                            const real_type theta4 , const real_type theta5 )
        {
          const int qubits[2] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
          const RotationTFXY< T > tmp( theta0 , theta1 , theta2 ,
                                       theta3 , theta4 , theta5 ) ;
          update( tmp.a() , tmp.b() , tmp.c() , tmp.d() ) ;
        } // RotationTFXYMatrix(qubit0,qubit1,theta0,theta1,theta2,theta3,theta4,theta5)

        /**
         * \brief Constructs a TFXY-rotation gate from the given TFXY-rotation
         *        gate `gate`.
         */
        RotationTFXYMatrix( const RotationTFXY< T >& gate )
        {
          setQubits( &(gate.qubits()[0]) ) ;
          update( gate.a() , gate.b() , gate.c() , gate.d() ) ;
        } // RotationTFXYMatrix(gate)

        /**
         * \brief Constructs a TFXY-rotation gate from the given TFXZ-rotation
         *        gate `gate`.
         */
        RotationTFXYMatrix( const RotationTFXZ< T >& gate )
        {
          setQubits( &(gate.qubits()[0]) ) ;
          update( gate.a() , gate.b() , gate.c() , gate.d() ) ;
        } // RotationTFXYMatrix(gate)

        /**
         * \brief Constructs a TFXY-rotation gate from the given TFYZ-rotation
         *        gate `gate`.
         */
        RotationTFXYMatrix( const RotationTFYZ< T >& gate )
        {
          setQubits( &(gate.qubits()[0]) ) ;
          update( gate.a() , gate.b() , gate.c() , gate.d() ) ;
        } // RotationTFXYMatrix(gate)

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
          using M = qclab::dense::SquareMatrix< T > ;
          return M( v_[0] ,   0   ,          0          , -std::conj( v_[3] ) ,
                      0   , v_[1] , -std::conj( v_[2] ) ,          0          ,
                      0   , v_[2] ,  std::conj( v_[1] ) ,          0          ,
                    v_[3] ,   0   ,          0          ,  std::conj( v_[0] ) );
        }

        // apply

        // print

        /**
         * \brief Writes the QASM code of this TFXY-rotation gate to the given
         *        `stream`.
         */
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          f3c::qgates::RotationTFXY< T >  tmp( *this ) ;
          return tmp.toQASM( stream , offset ) ;
        }

        // operator==

        // operator!=

        /// Checks if `other` equals this TFXY-rotation gate.
        inline bool equals( const qclab::QObject< T >& other ) const override {
          using TFXY = RotationTFXYMatrix< T > ;
          if ( const TFXY* p = dynamic_cast< const TFXY* >( &other ) ) {
            return ( p->a() == this->a() ) && ( p->b() == this->b() ) &&
                   ( p->c() == this->c() ) && ( p->d() == this->d() ) ;
          }
          return false ;
        }

        /// Returns the numerical value \f$a\f$ of this TFXY-rotation gate.
        inline T a() const { return v_[0] ; }

        /// Returns the numerical value \f$b\f$ of this TFXY-rotation gate.
        inline T b() const { return v_[1] ; }

        /// Returns the numerical value \f$c\f$ of this TFXY-rotation gate.
        inline T c() const { return v_[2] ; }

        /// Returns the numerical value \f$d\f$ of this TFXY-rotation gate.
        inline T d() const { return v_[3] ; }

        /// Returns the numerical values of this TFXY-rotation gate.
        inline const std::array< T , 4 >& values() const { return v_ ; }

        /// Updates this TFXY-rotation gate with the given parameters.
        inline void update( const T a , const T b , const T c , const T d ) {
          //const real_type eps = std::numeric_limits< real_type >::epsilon() ;
          //assert( std::abs( std::abs(a) * std::abs(a) +
          //                  std::abs(d) * std::abs(d) - T(1) ) < 10*eps ) ;
          //assert( std::abs( std::abs(b) * std::abs(b) +
          //                  std::abs(c) * std::abs(c) - T(1) ) < 10*eps ) ;
          v_[0] = a ;
          v_[1] = b ;
          v_[2] = c ;
          v_[3] = d ;
        }

        /// Multiplies `rhs` to this 2-qubit TFXY-rotation gate.
        inline RotationTFXYMatrix< T >& operator*=(
                                          const RotationTFXYMatrix< T >& rhs ) {
          assert( this->qubits()[0] == rhs.qubits()[0] ) ;
          assert( this->qubits()[1] == rhs.qubits()[1] ) ;
          update( rhs.a() * a() - std::conj( rhs.d() ) * d() ,
                  rhs.b() * b() - std::conj( rhs.c() ) * c() ,
                  rhs.c() * b() + std::conj( rhs.b() ) * c() ,
                  rhs.d() * a() + std::conj( rhs.a() ) * d() ) ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationTFXYMatrix< T > operator*( RotationTFXYMatrix< T > lhs ,
                                          const RotationTFXYMatrix< T >& rhs ) {
          assert( lhs.qubits()[0] == rhs.qubits()[0] ) ;
          assert( lhs.qubits()[1] == rhs.qubits()[1] ) ;
          lhs *= rhs ;
          return lhs ;
        }

      protected:
        /// Qubits of this TFXY-rotation gate.
        std::array< int , 2 >  qubits_ ;
        /// Complex values of this TFXY-rotation gate.
        std::array< T , 4 >    v_ ;

    } ; // class RotationTFXYMatrix

  } // namespace qgates

} // namespace f3c

#endif

