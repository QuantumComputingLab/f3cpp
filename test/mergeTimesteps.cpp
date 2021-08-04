#include "qclab/QCircuit.hpp"
#include "f3c/SquareCircuit.hpp"
#include "f3c/TriangleCircuit.hpp"
#include "f3c/qgates/functors.hpp"
#include <random>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

template <typename F>
auto triangleCircuit( const int n ) {

  using T = typename F::value_type ;
  using G = typename F::gate_type ;

  // random angles
  std::random_device                rd ;
  std::mt19937                      gen( rd() ) ;
  std::uniform_real_distribution<>  dis( 0.0, 1.0 ) ;

  // construct layers
  f3c::TriangleCircuit< T , G >  circuit( n ) ;
  int c = 0 ;
  for ( int l = 0; l < n-1; l++ ) {
    for ( int i = 0; i < n-l-1; i++ ) {
      circuit[c] = F::template init( n-i-2 , dis , gen ) ;
      c++ ;
    }
  }

  // return triangle circuit
  return circuit ;

}


template <typename F>
void timestep( qclab::QCircuit< typename F::value_type , typename F::gate_type >& layer1 ,
               qclab::QCircuit< typename F::value_type , typename F::gate_type >& layer2 ) {

  using G = typename F::gate_type ;

  // random angles
  std::random_device                rd ;
  std::mt19937                      gen( rd() ) ;
  std::uniform_real_distribution<>  dis( 0.0, 1.0 ) ;

  // update layers
  const int n = layer1.nbQubits() ;
  for ( int q = 0; q < n-1; q++ ) {
    if ( q % 2 == 0 ) {
      // even
      layer1[q/2] = F::template init( q , dis , gen ) ;
    } else {
      // odd
      layer2[q/2] = F::template init( q , dis , gen ) ;
    }
  }

}


template <typename F>
int timings( const int N , std::vector< int >& n , const int repeat ) {

  using G = typename F::gate_type ;
  using T = typename F::value_type ;
  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;
  const R tol = 100*eps ;

  // problem dimension
  std::cout << "qubits = " << N ;
#ifdef _OPENMP
  std::cout << ", omp_get_max_threads() = " << omp_get_max_threads() ;
#endif
  std::cout << std::endl << std::endl ;

  // random timestep
  qclab::QCircuit< T , G >  layer1( N , 0 , N/2 ) ;
  qclab::QCircuit< T , G >  layer2( N , 0 , (N-1)/2 ) ;

  // max norms
  std::vector< double >  nrms( n.size() , 0 ) ;

  // repeat
  for ( int r = 0; r < repeat; r++ ) {

    std::cout << "repeat = " << r+1 << ":" << std::endl ;

    // random triangle circuit
    auto triangle = triangleCircuit< F >( N ) ;

    // control circuit = triangle * timesteps
    qclab::QCircuit< T , G >  circuit( N ) ;
    for ( auto it = triangle.begin(); it != triangle.end(); ++it ) {
      circuit.push_back( std::make_unique< G >( **it ) ) ;
    }

    // merge timesteps
    int c = 0 ;
    for ( int i = 0; i < n.back(); i++ ) {

      // new random timestep
      timestep< F >( layer1 , layer2 ) ;

      // control circuit = triangle * timestep
      for ( auto it = layer1.begin(); it != layer1.end(); ++it ) {
        circuit.push_back( std::make_unique< G >( **it ) ) ;
      }
      for ( auto it = layer2.begin(); it != layer2.end(); ++it ) {
        circuit.push_back( std::make_unique< G >( **it ) ) ;
      }

      // merge timestep
      #pragma omp parallel for
      for ( size_t j = 0; j < layer1.nbGates(); j++ ) {
        triangle.merge( qclab::Side::Right , layer1[j] ) ;
      }
      #pragma omp parallel for
      for ( size_t j = 0; j < layer2.nbGates(); j++ ) {
        triangle.merge( qclab::Side::Right , layer2[j] ) ;
      }

      // Frobenius norm
      if ( n[c] == i+1 ) {
        std::cout << "timestep = " << i+1 << ":" << std::endl ;
        f3c::TriangleCircuit< T , G >  tmptriangle( triangle ) ;
        auto square = tmptriangle.toSquare() ;
        const double nrmF = qclab::nrmF( square , circuit ) ;
        if ( nrmF > nrms[c] ) nrms[c] = nrmF ;
        std::printf( "  * nrmF = %.4e  ( square: %3i , circuit: %5i )\n" ,
                     nrmF , int(square.nbGates()) , int(circuit.nbGates()) ) ;
        c++ ;
      }

    }

    std::cout << std::endl ;

  }

  // output
  std::cout << "results = [" << std::endl ;
  for ( size_t i = 0; i < nrms.size(); i++ ) {
    std::printf( "%6i, %10.4e" , n[i] , nrms[i] ) ;
    if ( i == nrms.size() - 1 ) {
      std::printf( "];\n" ) ;
    } else {
      std::printf( " ;\n" ) ;
    }
  }
  std::cout << std::endl ;

  // successful
  return 0 ;

}


int main( int argc , char *argv[] ) {

  int gate = 0 ;
  if ( argc <= 1 ) gate = 1 ;
  else if ( std::strcmp( argv[1] , "XY"   ) == 0 ) gate = 1 ;
  else if ( std::strcmp( argv[1] , "TFXY" ) == 0 ) gate = 2 ;

  char type   = 'd' ; if ( argc > 2 ) type   = argv[2][0] ;
  int  N      = 3 ;   if ( argc > 3 ) N      = std::stoi( argv[3] ) ;
  int  nmin   = 10 ;  if ( argc > 4 ) nmin   = std::stoi( argv[4] ) ;
  int  nmax   = 100 ; if ( argc > 5 ) nmax   = std::stoi( argv[5] ) ;
  int  step   = 10 ;  if ( argc > 6 ) step   = std::stoi( argv[6] ) ;
  int  repeat = 100 ; if ( argc > 7 ) repeat = std::stoi( argv[7] ) ;

  std::vector< int >  n ;
  if ( step > 0 ) {
    // linear scale
    for ( int i = nmin; i <= nmax; i += step ) {
      n.push_back( i ) ;
    }
    if ( n.back() != nmax ) { n.push_back( nmax ) ; }
  } else {
    // logarithmic scale
    std::vector< int > tmp( {
         1 ,    2 ,    3 ,    4 ,    5 ,    6 ,    7 ,    8 ,    9 ,    10 ,
        13 ,   16 ,   20 ,   25 ,   32 ,   40 ,   50 ,   63 ,   79 ,   100 ,
       126 ,  158 ,  200 ,  251 ,  316 ,  398 ,  501 ,  631 ,  794 ,  1000 ,
      1259 , 1585 , 1995 , 2512 , 3162 , 3981 , 5012 , 6310 , 7943 , 10000 } ) ;
    for ( int i = 0; i < tmp.size(); i++ ) {
      if ( nmin <= tmp[i] && tmp[i] <= nmax ) { n.push_back( tmp[i] ) ; }
    }
  }

  int r = 0 ;
  if ( gate == 1 ) {
    if ( type == 's' ) {
      // float
      std::cout << "\n+++ XY: merge timesteps (float) +++\n\n" ;
      using F = f3c::qgates::XYfunctor< std::complex< float > > ;
      r = timings< F >( N , n , repeat ) ;
      if ( r != 0 ) return r ;
    } else if ( type == 'd' ) {
      // double
      std::cout << "\n+++ XY: merge timesteps (double) +++\n\n" ;
      using F = f3c::qgates::XYfunctor< std::complex< double > > ;
      r = timings< F >( N , n , repeat ) ;
      if ( r != 0 ) return r + 1000 ;
    } else {
      return -2 ;
    }
  } else if ( gate == 2 ) {
    if ( type == 's' ) {
      // float
      std::cout << "\n+++ TFXY: merge timesteps (float) +++\n\n" ;
      using F = f3c::qgates::TFXYfunctor< std::complex< float > > ;
      r = timings< F >( N , n , repeat ) ;
      if ( r != 0 ) return r ;
    } else if ( type == 'd' ) {
      // double
      std::cout << "\n+++ TFXY: merge timesteps (double) +++\n\n" ;
      using F = f3c::qgates::TFXYfunctor< std::complex< double > > ;
      r = timings< F >( N , n , repeat ) ;
      if ( r != 0 ) return r + 1000 ;
    } else {
      return -3 ;
    }
  } else {
    return -1 ;
  }

  // successful
  std::cout << ">> end merge timesteps <<" << std::endl ;

  return 0 ;

}

