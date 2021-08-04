//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef f3c_TriangleCircuit_hpp
#define f3c_TriangleCircuit_hpp

#include "qclab/QCircuit.hpp"
#include "qclab/qgates/QGate2.hpp"
#include "f3c/turnover.hpp"

namespace f3c {

  template <typename T, typename G>
  class SquareCircuit ;

  /**
   * \class TriangleCircuit
   * \brief Class for representing a triangle quantum circuit.
   *
   * The quantum gates appear in the following pattern and ordering:
   *
   *       |
   *      | |
   *     | | |
   *    | | | |
   *   | | | | |
   *
   * We support 2 possible orderings of the gates
   *
   *   1) Ascending ordering:
   *
   *              04                  /
   *            03  08               / /
   *          02  07  11            / / /
   *        01  06  10  13         / / / /
   *      00  05  09  12  14      / / / / /
   *
   *
   *   2) Descending ordering:
   *
   *              10                  \
   *            06  11               \ \
   *          03  07  12            \ \ \
   *        01  04  08  13         \ \ \ \
   *      00  02  05  09  14      \ \ \ \ \
   *
   */
  template <typename T, typename G = qclab::qgates::QGate2< T >>
  class TriangleCircuit : public qclab::QCircuit< T , G >
  {

    public:
      /// Vector type of this triangle quantum circuit.
      using vector_type = typename qclab::QCircuit< T , G >::vector_type ;
      /// Size type of this triangle quantum circuit.
      using size_type   = typename qclab::QCircuit< T , G >::size_type ;

      /// Constructs a triangle quantum circuit of `nbQubits`.
      TriangleCircuit( const int nbQubits , const bool ascend = true )
      : qclab::QCircuit< T , G >( nbQubits , 0 , (nbQubits * (nbQubits-1)) / 2)
      , ascend_( ascend )
      {
        assert( nbQubits >= 2 ) ;
      } // TriangleCircuit(nbQubits,ascend)

      /// Copy constructor
      TriangleCircuit( const TriangleCircuit< T , G >& circuit )
      : qclab::QCircuit< T , G >( circuit.nbQubits() , circuit.offset() ,
                                  circuit.nbGates() )
      , ascend_( circuit.ascend() )
      {
        for ( size_t i = 0; i < circuit.nbGates(); i++ ) {
          auto gate = *circuit[i] ;
          this->gates_[i] = std::make_unique< decltype( gate ) >( gate ) ;
        }
      } // TriangleCircuit(circuit)

      /// Checks if this triangle quantum circuit is stored in ascending order.
      inline bool ascend() const { return ascend_ ; }

      /// Checks if this triangle quantum circuit is stored in descending order.
      inline bool descend() const { return !ascend_ ; }

      /// Reorders this triangle quantum circuit in ascending order.
      void makeAscend() {
        if ( ascend_ ) return ;
        // copy
        const auto nbGates = this->nbGates() ;
        vector_type gates( nbGates ) ;
        #pragma omp parallel for
        for ( size_type c = 0; c < nbGates; c++ ) {
          gates[c] = std::move( this->gates_[c] ) ;
        }
        // reodering
        const auto n = this->nbQubits() ;
        size_type c = 0 ;
        for ( int l = 0; l < n-1; l++ ) {
          for ( int i = 0; i < n-l-1; i++ ) {
            const auto idx = des2asc( l , i ) ;
            this->gates_[c] = std::move( gates[idx] ) ;
            c++ ;
          }
        }

        // update flag
        ascend_ = 1 ;
      }

      /// Reorders this triangle quantum circuit in descending order.
      void makeDescend() {
        if ( !ascend_ ) return ;
        // copy
        const auto nbGates = this->nbGates() ;
        vector_type gates( nbGates ) ;
        #pragma omp parallel for
        for ( size_type c = 0; c < nbGates; c++ ) {
          gates[c] = std::move( this->gates_[c] ) ;
        }
        // reodering
        const auto n = this->nbQubits() ;
        size_type c = 0 ;
        for ( int l = 0; l < n-1; l++ ) {
          for ( int i = 0; i < l+1; i++ ) {
            const auto idx = asc2des( l , i ) ;
            this->gates_[c] = std::move( gates[idx] ) ;
            c++ ;
          }
        }
        // update flag
        ascend_ = 0 ;
      }

      /**
       * \brief Merges the given gate `gate` on side `side` with this triangle
       *        quantum circuit.
       */
      void merge( qclab::Side side , std::unique_ptr< G >& gate ) {
        const int n = this->nbQubits() ;
        assert( gate->qubits()[0] < n - 1 ) ;
        assert( gate->qubits()[1] < n ) ;
        const int qubit = gate->qubit() ;
        auto& gates = this->gates_ ;
        std::unique_ptr< G >  gateA ;
        std::unique_ptr< G >  gateB ;
        std::unique_ptr< G >  gateC ;
        if ( ascend() ) {
          //
          // ascending
          //
          if ( side == qclab::Side::Left ) {
            // left
            int layer = 0 ;
            for ( int q = qubit; q < n-2 ; q++ ) {
              // turnovers
              const auto idx2 = ascIdx( layer , q + 1 ) ;
              const auto idx3 = ascIdx( layer , q     ) ;
              f3c::turnover( *gate , *gates[idx2] , *gates[idx3] ,
                             gateA , gateB , gateC ) ;
              gates[idx2] = std::move( gateA ) ;
              gates[idx3] = std::move( gateB ) ;
              gate = std::move( gateC ) ;
              layer++ ;
            }
            // fuse
            const auto idx = ascIdx( layer , n - 2 ) ;
            gates[idx] = std::make_unique< G >( (*gate) * (*gates[idx]) ) ;
          } else {
            // right
            const int layer = qubit ;
            for ( int q = qubit; q < n-2 ; q++ ) {
              // turnovers
              const auto idx1 = ascIdx( layer     , q     ) ;
              const auto idx2 = ascIdx( layer + 1 , q + 1 ) ;
              f3c::turnover( *gates[idx1] , *gates[idx2] , *gate ,
                             gateA , gateB , gateC ) ;
              gate = std::move( gateA ) ;
              gates[idx1] = std::move( gateB ) ;
              gates[idx2] = std::move( gateC ) ;
            }
            // fuse
            const auto idx = ascIdx( layer , n - 2 ) ;
            *gates[idx] *= *gate ;
          }
        } else {
          //
          // descending
          //
          if ( side == qclab::Side::Left ) {
            // left
            const int layer = n - qubit - 2 ;
            for ( int q = qubit; q < n-2 ; q++ ) {
              // turnovers
              const auto idx2 = desIdx( layer - 1 , q + 1 ) ;
              const auto idx3 = desIdx( layer     , q     ) ;
              f3c::turnover( *gate , *gates[idx2] , *gates[idx3] ,
                             gateA , gateB , gateC ) ;
              gates[idx2] = std::move( gateA ) ;
              gates[idx3] = std::move( gateB ) ;
              gate = std::move( gateC ) ;
            }
            // fuse
            const auto idx = desIdx( layer , n - 2 ) ;
            gates[idx] = std::make_unique< G >( (*gate) * (*gates[idx]) ) ;
          } else {
            // right
            int layer = n - 2 ;
            for ( int q = qubit; q < n-2 ; q++ ) {
              // turnovers
              const auto idx1 = desIdx( layer , q     ) ;
              const auto idx2 = desIdx( layer , q + 1 ) ;
              f3c::turnover( *gates[idx1] , *gates[idx2] , *gate ,
                             gateA , gateB , gateC ) ;
              gate = std::move( gateA ) ;
              gates[idx1] = std::move( gateB ) ;
              gates[idx2] = std::move( gateC ) ;
              layer-- ;
            }
            // fuse
            const auto idx = desIdx( layer , n - 2 ) ;
            *gates[idx] *= *gate ;
          }
        }
      }

      /**
       * \brief Merges the given gate `gate` on side `side` with this triangle
       *        quantum circuit.
       */
      void merge( qclab::Side side , const G& gate ) {
        std::unique_ptr< G >  p = std::make_unique< G >( gate ) ;
        merge( side , p ) ;
      }

      /// Converts this triangle quantum circuit into a square quantum circuit.
      SquareCircuit< T , G > toSquare() {
        makeAscend() ;
        const auto n = this->nbQubits() ;
        SquareCircuit< T , G >  square( n ) ;
        auto& gates = this->gates_ ;
        if ( n % 2 == 0 ) {
          //
          // even
          //
          const size_type stride = n/2 - 1 ;
          // turnovers
          for ( int i = 0; i < n/2 - 1; i++ ) {
            const auto l = n - 2*i - 2 ;
            const auto first = ascIdx( l , n - 2 ) ;
            for ( int j = 0; j <= 2*i; j++ ) {
              turnovers( l , n - 2 - j , gates[ first + j ] ) ;
              square[ i + j * stride ] = std::move( gates[ first + j ] ) ;
            }
          }
          // copy diagonal
          const size_type last = square.lastIdx( 0 ) ;
          #pragma omp parallel for
          for ( size_type i = 0; i < n-1; i++ ) {
            square[ last + i * stride ] = std::move( gates[i] ) ;
          }
          // copy subdiagonals
          #pragma omp parallel for
          for ( int l = 1; l < n-1; l += 2 ) {
            const auto idx = ascIdx( l , n - 2 ) ;
            const auto last = square.lastIdx( l + 1 ) ;
            for ( size_type i = 0; i < n-l-1; i++ ) {
              square[ last + i * stride ] = std::move( gates[ idx + i ] ) ;
            }
          }
        } else {
          //
          // odd
          //
          const size_type stride = n/2 ;
          // turnovers
          for ( int i = 0; i < n/2; i++ ) {
            const auto l = n - 2*i - 2 ;
            const auto first = ascIdx( l , n - 2 ) ;
            for ( int j = 0; j <= 2*i; j++ ) {
              turnovers( l , n - 2 - j , gates[ first + j ] ) ;
              const auto offset = j/2 * stride + ( j + 1 )/2 * ( stride - 1 ) ;
              square[ i + offset ] = std::move( gates[ first + j ] ) ;
            }
          }
          // copy (sub)diagonals
          #pragma omp parallel for
          for ( int l = 0; l < n-1; l += 2 ) {
            const auto idx = ascIdx( l , n - 2 ) ;
            const auto last = square.lastIdx( l + 1 ) ;
            for ( size_type i = 0; i < n-l-1; i++ ) {
              const auto offset = ( i + 1 )/2 * stride + i/2 * ( stride - 1 ) ;
              square[ last + offset ] = std::move( gates[ idx + i ] ) ;
            }
          }
        }
        return square ;
      }

      /// Returns the ascending to descending index.
      inline size_type asc2des( const int layer , const int index ) const {
        const auto n = this->nbQubits() ;
        assert( 0 <= layer ) ; assert( layer <= n - 2 ) ;
        assert( 0 <= index ) ; assert( index <= layer ) ;
        return layer + ( ( n - 1 ) * ( n - 2 ) ) / 2
                     - ( ( n - 1 - index ) * ( n - 2 - index ) ) / 2 ;
      }

      /// Returns the descending to ascending index.
      inline size_type des2asc( const int layer , const int index ) const {
        assert( 0 <= layer ) ; assert( layer <= this->nbQubits() - 2 ) ;
        assert( 0 <= index ) ; assert( index <= this->nbQubits() - 2 - layer ) ;
        return ( ( layer + 1 ) * ( layer + 2 ) ) / 2 - 1
                     + ( index * ( index + 1 ) ) / 2 + layer * index ;
      }

      /// Returns the ascending linear index.
      inline size_type ascIdx( const int layer , const int qubit ) const {
        assert( ascend() ) ;
        const auto n = this->nbQubits() ;
        assert( 0 <= layer ) ; assert( layer <= n - 2 ) ;
        assert( layer <= qubit ) ; assert( qubit <= n - 2 ) ;
        return ( layer * ( 2*n - layer - 1 ) ) / 2 + n - qubit - 2 ;
      }

      /// Returns the descending linear index.
      inline size_type desIdx( const int layer , const int qubit ) const {
        assert( descend() ) ;
        const auto n = this->nbQubits() ;
        assert( 0 <= layer ) ; assert( layer <= n - 2 ) ;
        assert( n - layer - 2 <= qubit ) ; assert( qubit <= n - 2 ) ;
        return ( layer * ( layer + 1 ) ) / 2 + layer - n + qubit + 2 ;
      }

    private:
      inline void turnovers( const int l , const int q ,
                             std::unique_ptr< G >& gate3 ) {
        std::unique_ptr< G >  gateA ;
        std::unique_ptr< G >  gateB ;
        std::unique_ptr< G >  gateC ;
        auto& gates = this->gates_ ;
        for ( int k = 0; k < l; k++ ) {
          const auto idx = ascIdx( l - k - 1 , q - k ) ;
          f3c::turnover( *gates[idx] , *gates[idx+1] , *gate3 ,
                         gateA , gateB , gateC ) ;
          gate3 = std::move( gateA ) ;
          gates[ idx     ] = std::move( gateB ) ;
          gates[ idx + 1 ] = std::move( gateC ) ;
        }
      }

    protected:
      bool  ascend_ ;  ///< Ordering of this triangle quantum circuit.

  } ; // class TriangleCircuit

} // namespace f3c

#endif

