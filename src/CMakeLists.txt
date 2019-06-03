########################
# Src and header files #
########################
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/include/nornir/)
include_directories(${CMAKE_SOURCE_DIR}/include/nornir/external/fastflow)
include_directories(${CMAKE_SOURCE_DIR}/src/external/mammut/)
include_directories(${CMAKE_SOURCE_DIR}/src/external/riff/src/)

file(GLOB SOURCES "*.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_SOURCE_DIR}/src/ompt.cpp")
 
####################
# Library building #
####################
add_library(nornir SHARED ${SOURCES})
#add_library(nornir_static STATIC ${SOURCES})

# Generate protocol fields names before building
add_custom_target(
    compile_ext
    COMMAND make -C ${CMAKE_SOURCE_DIR}/src/external
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
add_dependencies(nornir compile_ext)
#add_dependencies(nornir_static compile_ext)

install(TARGETS nornir 
        LIBRARY DESTINATION lib)

#install(TARGETS nornir_static 
#        ARCHIVE DESTINATION lib)

install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/nornir
        DESTINATION include)

####################
# Uninstall target #
####################
if(NOT TARGET uninstall)
    configure_file(
        "${CMAKE_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake)
endif()
