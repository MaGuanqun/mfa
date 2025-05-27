include_directories         (${CMAKE_SOURCE_DIR}/examples)          # needed for block.hpp and opts.h
include_directories         (${CMAKE_SOURCE_DIR}/examples/utility)

add_executable              (write_vtk              write_vtk.cpp writer.cpp)
add_executable              (write_two_value              write_two_value.cpp writer.cpp)
add_executable              (write_h              write_h.cpp writer.cpp)

target_link_libraries       (write_vtk              ${libraries})
target_link_libraries       (write_two_value              ${libraries})
target_link_libraries       (write_h              ${libraries})


set_property                (TARGET write_vtk  write_two_value write_h            PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

install                     (TARGETS write_vtk write_two_value write_h
                            DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/convert/
                            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_WRITE GROUP_EXECUTE
                            WORLD_READ WORLD_WRITE WORLD_EXECUTE)
