add_executable( f3c_time_evolution_XY timeEvolutionXY.cpp )
target_link_libraries( f3c_time_evolution_XY PUBLIC f3cpp qclabpp )

add_executable( f3c_time_evolution_XZ timeEvolutionXZ.cpp )
target_link_libraries( f3c_time_evolution_XZ PUBLIC f3cpp qclabpp )

add_executable( f3c_time_evolution_YZ timeEvolutionYZ.cpp )
target_link_libraries( f3c_time_evolution_YZ PUBLIC f3cpp qclabpp )

add_executable( f3c_time_evolution_TFXY timeEvolutionTFXY.cpp )
target_link_libraries( f3c_time_evolution_TFXY PUBLIC f3cpp qclabpp )

add_executable( f3c_time_evolution_TFXZ timeEvolutionTFXZ.cpp )
target_link_libraries( f3c_time_evolution_TFXZ PUBLIC f3cpp qclabpp )

add_executable( f3c_time_evolution_TFYZ timeEvolutionTFYZ.cpp )
target_link_libraries( f3c_time_evolution_TFYZ PUBLIC f3cpp qclabpp )
