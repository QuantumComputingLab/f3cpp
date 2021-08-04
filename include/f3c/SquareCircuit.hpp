//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_SquareCircuit_hpp
#define f3c_SquareCircuit_hpp

#include "qclab/QCircuit.hpp"
#include "qclab/qgates/QGate2.hpp"
#include "f3c/turnover.hpp"

namespace f3c {

  template <typename T, typename G>
  class TriangleCircuit ;

  /**
   * \class SquareCircuit
   * \brief Class for representing a square quantum circuit.
   *
   * The quantum gates appear in the following pattern and ordering:
   *
   *   A) Even number of qubits
   *
   *       |   |   |           00      05      10
   *         |   |   |             03      08      13
   *       |   |   |           01      06      11
   *         |   |   |             04      09      14
   *       |   |   |           02      07      12
   *
   *
   *   B) Odd number of qubits
   *
   *       |   |   |   |       00      06      12      18
   *         |   |   |             03      09      15
   *       |   |   |   |       01      07      13      19
   *         |   |   |             04      10      16
   *       |   |   |   |       02      08      14      20
   *         |   |   |             05      11      17
   *
   */
  template <typename T, typename G = qclab::qgates::QGate2< T >>
  class SquareCircuit : public qclab::QCircuit< T , G >
  {

    public:
      /// Size type of this square quantum circuit.
      using size_type = typename qclab::QCircuit< T , G >::size_type ;

      /// Constructs a square quantum circuit of `nbQubits`.
      SquareCircuit( const int nbQubits )
      : qclab::QCircuit< T , G >( nbQubits , 0 , (nbQubits * (nbQubits-1)) / 2)
      {
        assert( nbQubits >= 2 ) ;
      } // SquareCircuit(nbQubits)

      /// Copy constructor
      SquareCircuit( const SquareCircuit< T , G >& circuit )
      : qclab::QCircuit< T , G >( circuit.nbQubits() , circuit.offset() ,
                                  circuit.nbGates() )
      {
        for ( size_t i = 0; i < circuit.nbGates(); i++ ) {
          auto gate = *circuit[i] ;
          this->gates_[i] = std::make_unique< decltype( gate ) >( gate ) ;
        }
      } // SquareCircuit(circuit)

      /// Converts this square quantum circuit into a triangle quantum circuit.
      TriangleCircuit< T , G > toTriangle() {
        const auto n = this->nbQubits() ;
        TriangleCircuit< T , G >  triangle( n ) ;
        auto& gates = this->gates_ ;
        if ( n % 2 == 0 ) {
          //
          // even
          //
          const size_type stride = n/2 - 1 ;
          // copy diagonal
          const size_type last = lastIdx( 0 ) ;
          #pragma omp parallel for
          for ( size_type i = 0; i < n-1; i++ ) {
            triangle[i] = std::move( gates[ last + i * stride ] ) ;
          }
          // copy subdiagonals
          #pragma omp parallel for
          for ( int l = 1; l < n-1; l += 2 ) {
            const auto idx = triangle.ascIdx( l , n - 2 ) ;
            const auto last = lastIdx( l + 1 ) ;
            for ( size_type i = 0; i < n-l-1; i++ ) {
              triangle[ idx + i ] = std::move( gates[ last + i * stride ] ) ;
            }
          }
          // turnovers
          for ( int i = 0; i < n/2 - 1; i++ ) {
            const int sqLayer = n - 2*i - 4 ;
            const auto first = firstIdx( sqLayer ) ;
            for ( int q = 0; q <= sqLayer; q++ ) {
              turnovers( q , 2*i + 2 , gates[ first - q * stride ] , triangle );
              const auto idx = triangle.ascIdx( 2*i + 2 , 2*i + 2 + q ) ;
              triangle[idx] = std::move( gates[ first - q * stride ] ) ;
            }
          }
        } else {
          //
          // odd
          //
          const size_type stride = n/2 ;
          // copy (sub)diagonals
          #pragma omp parallel for
          for ( int l = 0; l < n-1; l += 2 ) {
            const auto idx = triangle.ascIdx( l , n - 2 ) ;
            const auto last = lastIdx( l + 1 ) ;
            for ( size_type i = 0; i < n-l-1; i++ ) {
              const auto offset = ( i + 1 )/2 * stride + i/2 * ( stride - 1 ) ;
              triangle[ idx + i ] = std::move( gates[ last + offset ] ) ;
            }
          }
          // turnovers
          for ( int i = 0; i < n/2; i++ ) {
            const int sqLayer = n - 2*i - 3 ;
            const auto first = firstIdx( sqLayer ) ;
            for ( int q = 0; q <= sqLayer; q++ ) {
              const auto offset = ( q + 1 )/2 * stride + q/2 * ( stride - 1 ) ;
              turnovers( q , 2*i + 1 , gates[ first - offset ] , triangle ) ;
              const auto idx = triangle.ascIdx( 2*i + 1 , 2*i + 1 + q ) ;
              triangle[idx] = std::move( gates[ first - offset ] ) ;
            }
          }
        }
        return triangle ;
      }

      /// Returns the linear index of the first gate in the given layer `layer`.
      inline size_type firstIdx( const size_type layer ) const {
        const auto n = this->nbQubits() ;
        assert( 0 <= layer ) ; assert( layer < n ) ;
        return layer/2 * ( n - 1 ) + ( layer % 2 ) * n/2 ;
      }

      /// Returns the linear index of the last gate in the given layer `layer`.
      inline size_type lastIdx( const size_type layer ) const {
        const auto n = this->nbQubits() ;
        assert( 0 <= layer ) ; assert( layer < n ) ;
        return n/2 - 1 + layer/2 * ( n - 1 ) + ( layer % 2 ) * ( (n+1)/2 - 1 ) ;
      }

    private:
      inline void turnovers( const int q , const int nb ,
                             std::unique_ptr< G >& gate1 ,
                             TriangleCircuit< T , G >& triangle ) {
        std::unique_ptr< G >  gateA ;
        std::unique_ptr< G >  gateB ;
        std::unique_ptr< G >  gateC ;
        for ( int l = 0; l < nb; l++ ) {
          const auto idx = triangle.ascIdx( l , q + l + 1 ) ;
          f3c::turnover( *gate1 , *triangle[ idx ] , *triangle[ idx + 1 ] ,
                         gateA , gateB , gateC ) ;
          triangle[ idx     ] = std::move( gateA ) ;
          triangle[ idx + 1 ] = std::move( gateB ) ;
          gate1 = std::move( gateC ) ;
        }
      }

  } ; // class SquareCircuit

} // namespace f3c

#endif

