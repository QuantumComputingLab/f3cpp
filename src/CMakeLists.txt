add_library( f3cpp io/INIFile.cpp )
target_include_directories( f3cpp PUBLIC ${PROJECT_SOURCE_DIR}/include )
target_compile_features( f3cpp PUBLIC cxx_std_17 )
if( TARGET OpenMP::OpenMP_CXX )
  target_link_libraries( f3cpp PUBLIC OpenMP::OpenMP_CXX )
endif()

