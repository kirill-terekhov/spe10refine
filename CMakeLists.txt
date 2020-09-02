cmake_minimum_required(VERSION 2.6)
project(SPE10)



find_package(inmost REQUIRED)

message(${INMOST_LIBRARY_DIRS})

link_directories(${INMOST_LIBRARY_DIRS})
include_directories(${INMOST_INCLUDE_DIRS})
add_definitions(${INMOST_DEFINITIONS})

add_executable(spe10grid spe10grid.cpp)

set_target_properties(spe10grid PROPERTIES COMPILE_FLAGS "${INMOST_COMPILE_FLAGS}")
set_target_properties(spe10grid PROPERTIES LINK_FLAGS "${INMOST_LINK_FLAGS}")
#target_include_directories(spe10grid PRIVATE $<TARGET_PROPERTY:AdaptiveMeshLib,INTERFACE_INCLUDE_DIRECTORIES>)

target_link_libraries(spe10grid ${INMOST_LIBRARIES} AdaptiveMeshLib)


install(TARGETS spe10grid EXPORT SPE10-targets RUNTIME DESTINATION bin/spe10)
install(FILES spe_perm.dat spe_phi.dat DESTINATION bin/spe10)

