#include "f3c/TriangleCircuit.hpp"
#include "f3c/qgates/functors.hpp"
#include <chrono>
#include <random>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

using TP = std::chrono::time_point< std::chrono::high_resolution_clock > ;

void tic( TP& t_bgn ) {
  t_bgn = std::chrono::high_resolution_clock::now() ;
} // tic

void tic( const std::string name , TP& t_bgn ) {
  std::cout << name << std::endl ;
  t_bgn = std::chrono::high_resolution_clock::now() ;
} // tic

double toc( const TP& t_bgn ) {
  TP t_end = std::chrono::high_resolution_clock::now() ;
  double time = std::chrono::duration< double >(t_end - t_bgn).count() ;
  return time ;
} // toc

double toc( const TP& t_bgn , TP& t_end ) {
  t_end = std::chrono::high_resolution_clock::now() ;
  double time = std::chrono::duration< double >(t_end - t_bgn).count() ;
  std::cout << "                                   " << time << "s\n" ;
  return time ;
} // toc


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
int timings( const std::vector< int >& N , const int outer , const int inner ,
             const int check ) {

  using G = typename F::gate_type ;
  using T = typename F::value_type ;
  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;
  const R tol = 100*eps ;

  // problem dimension
  std::cout << "qubits = " << N.front() << ":" << N.back()
            << ", outer = " << outer << ", inner = " << inner ;
#ifdef _OPENMP
  std::cout << ", omp_get_max_threads() = " << omp_get_max_threads() ;
#endif
  std::cout << std::endl << std::endl ;

  // timing variables
  TP  t_bgn ;
  TP  t_end ;

  std::vector< int >     qubits ;
  std::vector< double >  time ;
  for ( int n : N ) {

    std::cout << "n = " << n << ":" << std::endl ;
    qubits.push_back( n ) ;

    // random triangle circuit
    tic( "  * Constructing TriangleCircuit..." , t_bgn ) ;
    auto triangle = triangleCircuit< F >( n ) ;
    toc( t_bgn , t_end ) ;

    if ( check > 0 ) {

      // random timestep
      tic( "  * Constructing timestep..." , t_bgn ) ;
      qclab::QCircuit< T , G >  layer1( n , 0 , n/2 ) ;
      qclab::QCircuit< T , G >  layer2( n , 0 , (n-1)/2 ) ;
      timestep< F >( layer1 , layer2 ) ;
      toc( t_bgn , t_end ) ;

      // repeat
      for ( int i = 0; i < outer*inner; i++ ) {

        // control circuit = triangle * timestep
        qclab::QCircuit< T , G >  circuit( n ) ;
        if ( n <= check ) {
          tic( "  * Constructing ContoleCircuit..." , t_bgn ) ;
          for ( auto it = triangle.begin(); it != triangle.end(); ++it ) {
            circuit.push_back( std::make_unique< G >( **it ) ) ;
          }
          for ( auto it = layer1.begin(); it != layer1.end(); ++it ) {
            circuit.push_back( std::make_unique< G >( **it ) ) ;
          }
          for ( auto it = layer2.begin(); it != layer2.end(); ++it ) {
            circuit.push_back( std::make_unique< G >( **it ) ) ;
          }
          toc( t_bgn , t_end ) ;
        }

        // merge timestep
        tic( "  * Merge timestep..." , t_bgn ) ;
        #pragma omp parallel for
        for ( size_t i = 0; i < layer1.nbGates(); i++ ) {
          triangle.merge( qclab::Side::Right , layer1[i] ) ;
        }
        #pragma omp parallel for
        for ( size_t i = 0; i < layer2.nbGates(); i++ ) {
          triangle.merge( qclab::Side::Right , layer2[i] ) ;
        }
        const double ttot = toc( t_bgn , t_end ) ;
        if ( i == 0 ) time.push_back( ttot ) ;
        else if ( ttot < time.back() ) time.back() = ttot ;

        // check
        if ( n <= check ) {
          tic( "  * Frobenius norm..." , t_bgn ) ;
          const double nrmF = qclab::nrmF( triangle , circuit ) ;
          toc( t_bgn , t_end ) ;
          std::cout << "  * nrmF = " << nrmF << std::endl ;
          if ( nrmF > ( 1 << n ) * tol ) return n ;
        }
        std::cout << std::endl ;

        // new timestep
        if ( i < outer*inner - 1 ) {
          tic( "  * New timestep..." , t_bgn ) ;
          timestep< F >( layer1 , layer2 ) ;
          toc( t_bgn , t_end ) ;
        }

      }

    } else {

      // outer loop
      for ( int o = 0; o < outer; o++ ) {

        // random timesteps
        tic( "  * Constructing timesteps..." , t_bgn ) ;
        std::vector< qclab::QCircuit< T , G > >  layers1 ;
        std::vector< qclab::QCircuit< T , G > >  layers2 ;
        for ( int i = 0; i < inner; i++ ) {
          layers1.push_back( qclab::QCircuit< T , G >( n , 0 , n/2 ) ) ;
          layers2.push_back( qclab::QCircuit< T , G >( n , 0 , (n-1)/2 ) ) ;
          timestep< F >( layers1[i] , layers2[i] ) ;
        }
        const size_t nb1 = layers1[0].nbGates() ;
        const size_t nb2 = layers2[0].nbGates() ;
        toc( t_bgn , t_end ) ;

        // merge timesteps
        tic( "  * Merge timesteps..." , t_bgn ) ;
        for ( int i = 0; i < inner; i++ ) {
          // inner loop
          #pragma omp parallel for
          for ( size_t j = 0; j < nb1; j++ ) {
            triangle.merge( qclab::Side::Right , layers1[i][j] ) ;
          }
          #pragma omp parallel for
          for ( size_t j = 0; j < nb2; j++ ) {
            triangle.merge( qclab::Side::Right , layers2[i][j] ) ;
          }
        }
        const double ttot = toc( t_bgn , t_end ) / inner ;
        if ( o == 0 ) time.push_back( ttot ) ;
        else if ( ttot < time.back() ) time.back() = ttot ;

      }

    }

  }

  // output
  std::cout << "results = [" << std::endl ;
  for ( size_t i = 0; i < time.size(); i++ ) {
    std::printf( "%6i, %10.4e" , qubits[i] , time[i] ) ;
    if ( i == time.size() - 1 ) {
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

  char type  = 'd' ; if ( argc > 2 ) type  = argv[2][0] ;
  int  nmin  = 3 ;   if ( argc > 3 ) nmin  = std::stoi( argv[3] ) ;
  int  nmax  = 20 ;  if ( argc > 4 ) nmax  = std::stoi( argv[4] ) ;
  int  step  = 1 ;   if ( argc > 5 ) step  = std::stoi( argv[5] ) ;
  int  outer = 1 ;   if ( argc > 6 ) outer = std::stoi( argv[6] ) ;
  int  inner = 1 ;   if ( argc > 7 ) inner = std::stoi( argv[7] ) ;
  int  check = 8 ;   if ( argc > 8 ) check = std::stoi( argv[8] ) ;

  std::vector< int >  N ;
  if ( step > 0 ) {
    // linear scale
    for ( int i = nmin; i <= nmax; i += step ) {
      N.push_back( i ) ;
    }
    if ( N.back() != nmax ) { N.push_back( nmax ) ; }
  } else {
    // logarithmic scale
    std::vector< int > tmp( {
                2 ,    3 ,    4 ,    5 ,    6 ,    7 ,    8 ,    9 ,    10 ,
        13 ,   16 ,   20 ,   25 ,   32 ,   40 ,   50 ,   63 ,   79 ,   100 ,
       126 ,  158 ,  200 ,  251 ,  316 ,  398 ,  501 ,  631 ,  794 ,  1000 ,
      1259 , 1585 , 1995 , 2512 , 3162 , 3981 , 5012 , 6310 , 7943 , 10000 } ) ;
    for ( int i = 0; i < tmp.size(); i++ ) {
      if ( nmin <= tmp[i] && tmp[i] <= nmax ) { N.push_back( tmp[i] ) ; }
    }
  }

  int r = 0 ;
  if ( gate == 1 ) {
    if ( type == 's' ) {
      // float
      std::cout << "\n+++ XY: merge timestep (float) +++\n\n" ;
      using F = f3c::qgates::XYfunctor< std::complex< float > > ;
      r = timings< F >( N , outer , inner , check ) ;
      if ( r != 0 ) return r ;
    } else if ( type == 'd' ) {
      // double
      std::cout << "\n+++ XY: merge timestep (double) +++\n\n" ;
      using F = f3c::qgates::XYfunctor< std::complex< double > > ;
      r = timings< F >( N , outer , inner , check ) ;
      if ( r != 0 ) return r + 1000 ;
    } else {
      return -2 ;
    }
  } else if ( gate == 2 ) {
    if ( type == 's' ) {
      // float
      std::cout << "\n+++ TFXY: merge timestep (float) +++\n\n" ;
      using F = f3c::qgates::TFXYfunctor< std::complex< float > > ;
      r = timings< F >( N , outer , inner , check ) ;
      if ( r != 0 ) return r ;
    } else if ( type == 'd' ) {
      // double
      std::cout << "\n+++ TFXY: merge timestep (double) +++\n\n" ;
      using F = f3c::qgates::TFXYfunctor< std::complex< double > > ;
      r = timings< F >( N , outer , inner , check ) ;
      if ( r != 0 ) return r + 1000 ;
    } else {
      return -3 ;
    }
  } else {
    return -1 ;
  }

  // successful
  std::cout << ">> end merge timestep <<" << std::endl ;

  return 0 ;

}

