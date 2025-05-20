include_directories         (${CMAKE_SOURCE_DIR}/examples)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/examples/critical_point)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/examples/utility)          # needed for block.hpp and opts.h

add_executable              (derivative_control_point        derivative_control_point.cpp)
add_executable              (convert_cpt_to_vtk        convert_cpt_to_vtk.cpp ${CMAKE_SOURCE_DIR}/examples/convert/writer.cpp)
add_executable              (compute_critical_point        compute_critical_point.cpp)


add_executable              (unique_root_from_csv        unique_root_from_csv.cpp)
add_executable              (convert_root_to_vtk        convert_root_to_vtk.cpp ${CMAKE_SOURCE_DIR}/examples/convert/writer.cpp)
add_executable              (test_time        test_time.cpp)
add_executable              (root_error        root_error.cpp ${CMAKE_SOURCE_DIR}/examples/convert/writer.cpp)


target_link_libraries       (derivative_control_point        ${libraries})
target_link_libraries       (convert_cpt_to_vtk        ${libraries})
target_link_libraries       (compute_critical_point        ${libraries})
# target_link_libraries       (compute_critical_point_jacobi_set        ${libraries})
target_link_libraries       (convert_root_to_vtk        ${libraries})
target_link_libraries       (test_time        ${libraries})
target_link_libraries       (root_error        ${libraries})
target_link_libraries       (unique_root_from_csv        ${libraries})

# set_property                (TARGET  derivative_control_point convert_cpt_to_vtk  compute_critical_point PROPERTY convert_root_to_vtk INSTALL_RPATH_USE_LINK_PATH TRUE)
# install                     (TARGETS derivative_control_point convert_cpt_to_vtk compute_critical_point convert_root_to_vtk
set_property                (TARGET derivative_control_point convert_cpt_to_vtk unique_root_from_csv convert_root_to_vtk test_time root_error compute_critical_point PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
install                     (TARGETS derivative_control_point convert_cpt_to_vtk unique_root_from_csv convert_root_to_vtk test_time root_error compute_critical_point DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/critical_point/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
