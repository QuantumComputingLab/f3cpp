add_executable( f3c_time_evolution_XY timeEvolutionXY.cpp )
target_link_libraries( f3c_time_evolution_XY PUBLIC f3cpp qclabpp )

add_executable( f3c_time_evolution_TFXY timeEvolutionTFXY.cpp )
target_link_libraries( f3c_time_evolution_TFXY PUBLIC f3cpp qclabpp )
