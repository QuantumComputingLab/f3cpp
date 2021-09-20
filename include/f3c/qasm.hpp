//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_qasm_hpp
#define f3c_qasm_hpp

#include "qclab/qasm.hpp"

namespace f3c {

  /**
   * \brief Returns a qasm string for a Rotation XY gate on the given qubits
   *        `qubit0` and `qubit1`, and with angles `theta0` and `theta1`.
   */
  template <typename T>
  inline auto qasmRxy( const int qubit0 , const int qubit1 ,
                       const T theta0 , const T theta1 ) {
    const T pi2 = 2 * std::atan(1) ;
    std::stringstream stream ;
    stream << qclab::qasmRx( qubit0 , pi2 )
           << qclab::qasmRx( qubit1 , pi2 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , theta0 )
           << qclab::qasmRz( qubit1 , theta1 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , -pi2 )
           << qclab::qasmRx( qubit1 , -pi2 ) ;
    return stream.str() ;
  }

  /**
   * \brief Returns a qasm string for a Rotation XZ gate on the given qubits
   *        `qubit0` and `qubit1`, and with angles `theta0` and `theta1`.
   */
  template <typename T>
  inline auto qasmRxz( const int qubit0 , const int qubit1 ,
                       const T theta0 , const T theta1 ) {
    std::stringstream stream ;
    stream << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , theta0 )
           << qclab::qasmRz( qubit1 , theta1 )
           << qclab::qasmCX( qubit0 , qubit1 ) ;
    return stream.str() ;
  }

  /**
   * \brief Returns a qasm string for a Rotation YZ gate on the given qubits
   *        `qubit0` and `qubit1`, and with angles `theta0` and `theta1`.
   */
  template <typename T>
  inline auto qasmRyz( const int qubit0 , const int qubit1 ,
                       const T theta0 , const T theta1 ) {
    const T pi2 = 2 * std::atan(1) ;
    std::stringstream stream ;
    stream << qclab::qasmRz( qubit0 , pi2 )
           << qclab::qasmRz( qubit1 , pi2 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , theta0 )
           << qclab::qasmRz( qubit1 , theta1 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRz( qubit0 , -pi2 )
           << qclab::qasmRz( qubit1 , -pi2 ) ;
    return stream.str() ;
  }


  /**
   * \brief Returns a qasm string for a Rotation TFXY gate on the given qubits
   *        `qubit0` and `qubit1`, and with angles `theta0`, `theta1`, `theta2`,
   *        `theta3`, `theta4`, and `theta5`.
   */
  template <typename T>
  inline auto qasmTFRxy( const int qubit0 , const int qubit1 ,
                         const T theta0 , const T theta1 , const T theta2 ,
                         const T theta3 , const T theta4 , const T theta5 ) {
    const T pi2 = 2 * std::atan(1) ;
    std::stringstream stream ;
    stream << qclab::qasmRz( qubit0 , theta0 )
           << qclab::qasmRz( qubit1 , theta1 )
           << qclab::qasmRx( qubit0 , pi2 )
           << qclab::qasmRx( qubit1 , pi2 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , theta2 )
           << qclab::qasmRz( qubit1 , theta3 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , -pi2 )
           << qclab::qasmRx( qubit1 , -pi2 )
           << qclab::qasmRz( qubit0 , theta4 )
           << qclab::qasmRz( qubit1 , theta5 ) ;
    return stream.str() ;
  }

  /**
   * \brief Returns a qasm string for a Rotation TFXZ gate on the given qubits
   *        `qubit0` and `qubit1`, and with angles `theta0`, `theta1`, `theta2`,
   *        `theta3`, `theta4`, and `theta5`.
   */
  template <typename T>
  inline auto qasmTFRxz( const int qubit0 , const int qubit1 ,
                         const T theta0 , const T theta1 , const T theta2 ,
                         const T theta3 , const T theta4 , const T theta5 ) {
    std::stringstream stream ;
    stream << qclab::qasmRy( qubit0 , theta0 )
           << qclab::qasmRy( qubit1 , theta1 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , theta2 )
           << qclab::qasmRz( qubit1 , theta3 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRy( qubit0 , theta4 )
           << qclab::qasmRy( qubit1 , theta5 ) ;
    return stream.str() ;
  }

  /**
   * \brief Returns a qasm string for a Rotation TFYZ gate on the given qubits
   *        `qubit0` and `qubit1`, and with angles `theta0`, `theta1`, `theta2`,
   *        `theta3`, `theta4`, and `theta5`.
   */
  template <typename T>
  inline auto qasmTFRyz( const int qubit0 , const int qubit1 ,
                         const T theta0 , const T theta1 , const T theta2 ,
                         const T theta3 , const T theta4 , const T theta5 ) {
    const T pi2 = 2 * std::atan(1) ;
    std::stringstream stream ;
    stream << qclab::qasmRx( qubit0 , theta0 )
           << qclab::qasmRx( qubit1 , theta1 )
           << qclab::qasmRz( qubit0 , pi2 )
           << qclab::qasmRz( qubit1 , pi2 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRx( qubit0 , theta2 )
           << qclab::qasmRz( qubit1 , theta3 )
           << qclab::qasmCX( qubit0 , qubit1 )
           << qclab::qasmRz( qubit0 , -pi2 )
           << qclab::qasmRz( qubit1 , -pi2 )
           << qclab::qasmRx( qubit0 , theta4 )
           << qclab::qasmRx( qubit1 , theta5 ) ;
    return stream.str() ;
  }

} // namespace f3c

#endif

