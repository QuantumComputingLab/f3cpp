//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_util_hpp
#define f3c_util_hpp

#include "qclab/util.hpp"
#include <cmath>
#include <array>
#include <complex>

/// Fast Free Fermion Compiler namespace
namespace f3c {

  /// Returns the signum function of the real number `val`.
  template <typename T>
  int sign( const T val ) {
    return ( T(0) < val ) - ( val < T(0) ) ;
  }

  /// Returns the Frobenius norm of a 2 x 2 matrix.
  template <typename R>
  inline R norm22( const std::complex< R >* A ) {

    // scale
    const double scale1 = std::max( std::abs( A[0] ) , std::abs( A[1] ) ) ;
    const double scale2 = std::max( std::abs( A[2] ) , std::abs( A[3] ) ) ;
    double scale = std::max( scale1 , scale2 ) ;
    if ( scale == 0.0 ) { scale = 1.0 ; }

    // trace^2
    double trace = 0 ;
    for ( int i = 0; i < 4; i++ ) {
      const double real = std::real( A[i] ) / scale ;
      const double imag = std::imag( A[i] ) / scale ;
      trace += real * real + imag * imag ;
    }

    // norm
    return scale * std::sqrt( trace ) ;

  }

  /// Multiplies two 2 x 2 matrices.
  template <typename T>
  inline void gemm22( const T* A , const T* B , T* C ) {

    C[0] = A[0]*B[0] + A[2]*B[1] ;
    C[1] = A[1]*B[0] + A[3]*B[1] ;
    C[2] = A[0]*B[2] + A[2]*B[3] ;
    C[3] = A[1]*B[2] + A[3]*B[3] ;

  }

  /**
   * \brief Computes the eigenvalues of a 2 x 2 matrix pencil.
   *
   * The eigenvalues are \f$\lambda_1 = \alpha_1/\beta\f$, and
   *                     \f$\lambda_2 = \alpha_2/\beta\f$.
   */
  template <typename T>
  void eig22( const T* A , const T* B , T& a1 , T& a2 , T& b ) {

    const T AB11 = A[0] * B[3] - A[1] * B[2];
    const T AB12 = A[2] * B[3] - A[3] * B[2];
    const T AB21 = A[1] * B[0] - A[0] * B[1];
    const T AB22 = A[3] * B[0] - A[2] * B[1];

    const T t = ( AB11 + AB22 ) / T(2) ;
    const T d = std::sqrt( t * t + AB12 * AB21 - AB11 * AB22 ) ;

    a1 = t + d ;
    a2 = t - d ;
    b = B[0] * B[3] - B[2] * B[1] ;

  }

  /**
   * \brief Computes a real rotation matrix (cos and sin) to introduce a zero
   *        in a vector of length 2.
   *
   * Implementation based on Chapter 1 of:
   *   "Core-Chasing Algorithms for the Eigenvalue Problem". Jared L. Aurentz,
   *   Thomas Mach, Leonardo Robol, Raf Vandebril, and David S. Watkins (2018).
   */
  template <typename R>
  inline std::tuple< R , R > rotateToZero( const R x , const R y ) {
    R c ;
    R s ;
    if ( y == 0 ) {
      if ( x == 0 ) {
        c = 1 ;
      } else {
        c = std::abs(x) / x ;
      }
      s = 0 ;
    } else {
      if ( std::abs(x) >= std::abs(y) ) {
        const auto theta = sign( x ) ;
        const auto t = y / x ;
        const auto r = std::sqrt( 1 + std::abs(t) * std::abs(t) ) ;
        c = theta / r ;
        s = t * c ;
      } else {
        const auto theta = sign( y ) ;
        const auto t = x / y ;
        const auto r = std::sqrt( 1 + std::abs(t) * std::abs(t) ) ;
        s = theta / r ;
        c = t * s ;
      }
    }
    return { c , s } ;
  }

  /**
   * \brief Computes a complex rotation matrix (cos and sin) to introduce a zero
   *        in a vector of length 2.
   *
   * Implementation based on Chapter 1 of:
   *   "Core-Chasing Algorithms for the Eigenvalue Problem". Jared L. Aurentz,
   *   Thomas Mach, Leonardo Robol, Raf Vandebril, and David S. Watkins (2018).
   */
  template <typename R>
  inline std::tuple< std::complex< R > , std::complex< R > >
    rotateToZero( const std::complex< R >& x , const std::complex< R >& y ) {

    std::complex< R >  c ;
    std::complex< R >  s ;
    if ( std::real(y) == 0 && std::imag(y) == 0 ) {
      if ( std::real(x) == 0 && std::imag(x) == 0 ) {
        c = 1 ;
      } else {
        c = std::abs(x) / x ;
      }
      s = 0 ;
    } else {
      if ( std::abs(x) >= std::abs(y) ) {
        const auto theta = std::conj( x / std::abs(x) ) ;
        const auto t = y / x ;
        const auto r = std::sqrt( R(1) + std::abs(t) * std::abs(t) ) ;
        c = theta / r ;
        s = std::conj(t) * c ;
      } else {
        const auto theta = std::conj( y / std::abs(y) ) ;
        const auto t = x / y ;
        const auto r = std::sqrt( R(1) + std::abs(t) * std::abs(t) ) ;
        s = theta / r ;
        c = std::conj(t) * s ;
      }
    }
    return { c , s } ;

  }

  /**
   * \brief Computes a complex rotation matrix to introduce a zero in a vector
   *        of length 2, i.e., G * [x; y] = [r; 0].
   */
  template <typename R>
  inline void rotateToZeroL( const std::complex< R >& x ,
                             const std::complex< R >& y ,
                             std::complex< R >* G ) {

    const auto [ c , s ] = rotateToZero( x , y ) ;

    // G = [       c       s  ]
    //     [ -conj(s) conj(c) ]
    G[0] = c ;
    G[1] = -std::conj(s) ;
    G[2] = s ;
    G[3] =  std::conj(c) ;

  }

  /**
   * \brief Computes a complex rotation matrix to introduce a zero in a vector
   *        of length 2, i.e., [y, x] * G = [0, r].
   */
  template <typename R>
  inline void rotateToZeroR( const std::complex< R >& x ,
                             const std::complex< R >& y ,
                             std::complex< R >* G ) {

    const auto [ c , s ] = rotateToZero( x , y ) ;

    // G = [  conj(c) s ]
    //     [ -conj(s) c ]
    G[0] =  std::conj(c) ;
    G[1] = -std::conj(s) ;
    G[2] = s ;
    G[3] = c ;

  }

  /**
   * \brief Computes a diagonalization for a 2 x 2 diagonalizable matrix pencil.
   *
   * The function is typically used for diagonalizable pencils of the form:
   * \f[\begin{bmatrix} a   &  b\\ c & d \end{bmatrix} - \lambda
   * \begin{bmatrix} \bar{d} & -\bar{c} \\ -\bar{b} & \bar{a} \end{bmatrix}\f]
   */
  template <typename T>
  void diagonalize22( T* A , T* B , T* Q , T* Z ) {

    // parameters
    int maxit = 100 ;
    int stepit = 3 ;

    // Q = Z = eye(2)
    Q[0] = 1 ; Q[1] = 0 ; Q[2] = 0 ; Q[3] = 1 ;
    Z[0] = 1 ; Z[1] = 0 ; Z[2] = 0 ; Z[3] = 1 ;

    // work array
    std::array< T , 4 >  V_ ;  T* V = V_.data() ;
    std::array< T , 4 >  W_ ;  T* W = W_.data() ;

    // norms and eps
    const auto nrmA = norm22( A ) ;
    const auto nrmB = norm22( B ) ;
    const auto eps = std::numeric_limits< qclab::real_t< T > >::epsilon() ;
    const auto tolA = 12 * eps * nrmA ;
    const auto tolB = 12 * eps * nrmB ;

    // loop
    T a1 ;
    T a2 ;
    T b ;
    for ( int it = 0; it < maxit; it++ ) {

      // compute eigenvalues
      if ( (it+1) % stepit == 0 ) {
        // exceptional shift to try to get out of stagnation
        if ( A[3] != T(0) ) {
          a1 = A[3] ;
          a2 = 0 ;
          b  = B[3] ;
        } else {
          a1 = A[0] ;
          a2 = 0 ;
          b  = B[0] ;
        }
      } else {
        eig22( A , B , a1 , a2 , b ) ;
      }

      // W = b*A + a*B
      if ( std::abs( a1 ) > std::abs( a2 ) ) {
        W[0] = b * A[0] - a1 * B[0] ;
        W[1] = b * A[1] - a1 * B[1] ;
        W[2] = b * A[2] - a1 * B[2] ;
        W[3] = b * A[3] - a1 * B[3] ;
      } else {
        W[0] = b * A[0] - a2 * B[0] ;
        W[1] = b * A[1] - a2 * B[1] ;
        W[2] = b * A[2] - a2 * B[2] ;
        W[3] = b * A[3] - a2 * B[3] ;
      }

      // rotate to 0 (left)
      if ( std::abs(W[0]) + std::abs(W[1]) > std::abs(W[2]) + std::abs(W[3]) ) {
        rotateToZeroL( W[0] , W[1] , V ) ;
      } else {
        rotateToZeroL( W[2] , W[3] , V ) ;
      }

      // A = V*A
      W[0] = A[0] ; W[1] = A[1] ; W[2] = A[2] ; W[3] = A[3] ;
      gemm22( V , W , A ) ;
      // B = V*B
      W[0] = B[0] ; W[1] = B[1] ; W[2] = B[2] ; W[3] = B[3] ;
      gemm22( V , W , B ) ;
      // Q = V*Q
      W[0] = Q[0] ; W[1] = Q[1] ; W[2] = Q[2] ; W[3] = Q[3] ;
      gemm22( V , W , Q ) ;

      // rotate to 0 (right)
      if ( std::abs(A[1]) + std::abs(A[3]) > std::abs(B[1]) + std::abs(B[3]) ) {
        rotateToZeroR( A[3] , A[1] , V ) ;
      } else {
        rotateToZeroR( B[3] , B[1] , V ) ;
      }

      // A = A*V
      W[0] = A[0] ; W[1] = A[1] ; W[2] = A[2] ; W[3] = A[3] ;
      gemm22( W , V , A ) ;
      // B = B*V
      W[0] = B[0] ; W[1] = B[1] ; W[2] = B[2] ; W[3] = B[3] ;
      gemm22( W , V , B ) ;
      // Z = Z*V
      W[0] = Z[0] ; W[1] = Z[1] ; W[2] = Z[2] ; W[3] = Z[3] ;
      gemm22( W , V , Z ) ;

      // check
      if ( ( std::abs( A[1] ) + std::abs( A[2] ) < tolA ) &&
           ( std::abs( B[1] ) + std::abs( B[2] ) < tolB ) ) {
        break ;
      }
      if ( it == maxit-1 ) {
        std::cout << "WARNING: maximum number of iterations in diagonalize22 "
                  << "reached, results may be inaccurate." << std::endl ;
        std::printf( "  |A(1,2)| + |A(2,1)| = %.4e  (tolA = %.4e)\n" ,
                     std::abs( A[1] ) + std::abs( A[2] ) , tolA ) ;
        std::printf( "  |B(1,2)| + |B(2,1)| = %.4e  (tolB = %.4e)\n" ,
                     std::abs( B[1] ) + std::abs( B[2] ) , tolB ) ;
      }

    }

    // A = diag(diag(Q*A*Z))
    A[1] = 0 ;
    A[2] = 0 ;
    // B = diag(diag(Q*B*Z))
    B[1] = 0 ;
    B[2] = 0 ;

  }

} // namespace f3c

#endif

