add_executable( f3c_tests run.cpp
                          qgates/RotationXY.cpp
                          qgates/RotationXZ.cpp
                          qgates/RotationYZ.cpp
                          qgates/RotationTFXY.cpp
                          qgates/RotationTFXZ.cpp
                          qgates/RotationTFYZ.cpp
                          qgates/RotationTFXYMatrix.cpp
                          SquareCircuit.cpp
                          TriangleCircuit.cpp
                          util.cpp
                          turnoverSU2.cpp
                          turnover.cpp
                          concepts.cpp
              )
target_link_libraries( f3c_tests PUBLIC f3cpp qclabpp gtest )
target_include_directories( f3c_tests PUBLIC ${PROJECT_SOURCE_DIR}/test )

add_executable( f3c_time_merge_timestep mergeTimestep.cpp )
target_link_libraries( f3c_time_merge_timestep PUBLIC f3cpp qclabpp )
target_include_directories( f3c_time_merge_timestep PUBLIC ${PROJECT_SOURCE_DIR}/test )

add_executable( f3c_time_merge_timesteps mergeTimesteps.cpp )
target_link_libraries( f3c_time_merge_timesteps PUBLIC f3cpp qclabpp )
target_include_directories( f3c_time_merge_timesteps PUBLIC ${PROJECT_SOURCE_DIR}/test )

add_executable( f3c_time_square2triangle square2triangle.cpp )
target_link_libraries( f3c_time_square2triangle PUBLIC f3cpp qclabpp )
target_include_directories( f3c_time_square2triangle PUBLIC ${PROJECT_SOURCE_DIR}/test )

add_executable( f3c_time_triangle2square triangle2square.cpp )
target_link_libraries( f3c_time_triangle2square PUBLIC f3cpp qclabpp )
target_include_directories( f3c_time_triangle2square PUBLIC ${PROJECT_SOURCE_DIR}/test )

