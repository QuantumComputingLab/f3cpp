#include "timeEvolution.hpp"
#include "f3c/io/INIFile.hpp"
#include "f3c/io/param.hpp"

int main( int argc , char *argv[] ) {

  using R = double ;
  using T = std::complex< R > ;
  using F = f3c::qgates::TFYZfunctor< T > ;
  using P = f3c::Param< R > ;
  using CV = f3c::ConstValue< R > ;

  // arguments
  if ( argc < 2 ) {
    std::cout << "ERROR: filename is required!" << std::endl ;
    return -1 ;
  }
  std::string filename( argv[1] ) ;
  {
    std::ifstream stream( filename ) ;
    if ( !stream.good() ) {
      std::cout << "ERROR: \"" << filename << "\" does not exist!" << std::endl;
      return -2 ;
    }
  }
  f3c::io::INIFile file( filename ) ;

  std::cout << "********************\n"
            << "***  TFYZ MODEL  ***\n"
            << "********************\n" << std::endl ;

  // qubits
  int N ;
  if ( file.contains( "Qubits.number" ) ) {
    N = file.value< int >( "Qubits.number" ) ;
  } else {
    std::cout << "ERROR: number of qubits is required!" << std::endl ;
    return -3 ;
  }

  // Trotter
  int n ;
  if ( file.contains( "Trotter.steps" ) ) {
    n = file.value< int >( "Trotter.steps" ) ;
  } else {
    std::cout << "ERROR: number of Trotter steps is required!" << std::endl ;
    return -4 ;
  }
  double dt ;
  if ( file.contains( "Trotter.dt" ) ) {
    dt = file.value< double >( "Trotter.dt" ) ;
  } else {
    std::cout << "ERROR: timestep size is required!" << std::endl ;
    return -5 ;
  }

  // parameter 0
  std::unique_ptr< P >  P0_ = std::make_unique< CV >( 0 ) ;
  P* P0 = P0_.get() ;

  // parameter hx
  std::unique_ptr< P >  hx_ ;
  if ( f3c::io::param( N , n , file , "hx" , hx_ ) != 0 ) return -6 ;
  P* hx = hx_.get() ;

  // parameter Jy
  std::unique_ptr< P >  Jy_ ;
  if ( f3c::io::param( N , n , file , "Jy" , Jy_ ) != 0 ) return -7 ;
  P* Jy = Jy_.get() ;

  // parameter Jz
  std::unique_ptr< P >  Jz_ ;
  if ( f3c::io::param( N , n , file , "Jz" , Jz_ ) != 0 ) return -8 ;
  P* Jz = Jz_.get() ;

  // Output
  std::string  name = "out" ;
  if ( file.contains( "Output.name" ) ) {
    name = file.value< std::string >( "Output.name" ) ;
  }
  int imin = 1 ;
  if ( file.contains( "Output.imin" ) ) {
    imin = file.value< int >( "Output.imin" ) ;
  }
  int imax = n ;
  if ( file.contains( "Output.imax" ) ) {
    imax = file.value< int >( "Output.imax" ) ;
  }
  int step = 1 ;
  if ( file.contains( "Output.step" ) ) {
    step = file.value< int >( "Output.step" ) ;
  }

  // Debug
  int debug = 0 ;
  if ( file.contains( "Debug.value" ) ) {
    debug = file.value< int >( "Debug.value" ) ;
  }

  std::cout << "  N = " << N << ": timesteps = " << n << " , dt = " << dt
            << "\n\n"
            << "    hx = " << *hx << "\n"
            << "    Jy = " << *Jy << "\n"
            << "    Jz = " << *Jz << "\n\n" ;

  return timeEvolution< F , P >( N , n , dt , imin , imax , step ,
                                 hx , P0 , P0 , P0 , Jy , Jz , name , debug ) ;

}

