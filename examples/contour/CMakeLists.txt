include_directories         (${CMAKE_SOURCE_DIR}/examples)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/examples/contour)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/examples/utility)  

add_executable              (compute_critical_point_jacobi_set        compute_critical_point_jacobi_set.cpp)
add_executable              (isocontour        isocontour.cpp)
add_executable              (jacobi_set        jacobi_set.cpp)
add_executable              (compute_function_error compute_function_error.cpp)

target_link_libraries         (compute_critical_point_jacobi_set ${libraries})
target_link_libraries       (isocontour        ${libraries})
target_link_libraries       (jacobi_set        ${libraries})
target_link_libraries       (compute_function_error ${libraries})

set_property                (TARGET compute_critical_point_jacobi_set jacobi_set compute_function_error isocontour PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install                     (TARGETS compute_critical_point_jacobi_set jacobi_set compute_function_error isocontour
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/contour/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)