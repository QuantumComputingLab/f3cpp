cmake_minimum_required( VERSION 3.16 FATAL_ERROR )

# project
project( F3CPP VERSION 0.1.0 LANGUAGES CXX )

# openmp
find_package( OpenMP )

# fetch content
include( FetchContent )

# qclabpp
FetchContent_Declare( qclabpp
  GIT_REPOSITORY https://github.com/QuantumComputingLab/qclabpp.git
)
FetchContent_MakeAvailable( qclabpp )

# GTest
FetchContent_Declare( gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
)
FetchContent_MakeAvailable( gtest )

# f3c/src
add_subdirectory( src )

# f3c/test
add_subdirectory( test )

# f3c/examples
add_subdirectory( examples )

# f3c/python
add_subdirectory( python )

# documentation
find_package( Doxygen )
if( DOXYGEN_FOUND )
  configure_file( ${CMAKE_SOURCE_DIR}/doc/doxygen/doxygen.dox.in
    ${CMAKE_BINARY_DIR}/doxygen.dox @ONLY)
  add_custom_target( f3c_doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doxygen.dox
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating API documentation with doxygen" VERBATIM )
endif()

