########################
# Src and header files #
########################
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/nornir/external/fastflow)

file(GLOB SOURCES "*.cpp" "*.c")
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/src/ompt.cpp")
 
####################
# Library building #
####################

# Must be set before add_library
set(EXTERNAL_INSTALL_LOCATION ${PROJECT_BINARY_DIR}/external)
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_library(nornir SHARED ${SOURCES})
add_library(nornir_static STATIC ${SOURCES})

# Download dependencies
include(ExternalProject)

ExternalProject_Add(mammut_repo
    GIT_REPOSITORY https://github.com/danieledesensi/mammut
    GIT_TAG 02f9c7d
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
)

ExternalProject_Add(riff_repo
    GIT_REPOSITORY https://github.com/danieledesensi/riff
    GIT_TAG 96aaab1
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
)

# Generate protocol fields names before building
add_custom_target(
    compile_ext
    COMMAND make -C ${PROJECT_SOURCE_DIR}/src/external
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

# Includes
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/nornir
        DESTINATION include)

# TODO bin


# Check if mlpack and gsl are present
find_library(MLPACK_LIBRARY NAMES mlpack)
find_library(GSL_LIBRARY NAMES gsl)

# Dynamic Library
add_dependencies(nornir compile_ext mammut_repo riff_repo)
target_link_libraries(nornir riff mammut ${MLPACK_LIBRARY} ${GSL_LIBRARY} pthread) # TODO: Check if dependencies are installed
install(TARGETS nornir LIBRARY DESTINATION lib)

# Static Library
add_dependencies(nornir_static compile_ext mammut_repo riff_repo)
target_link_libraries(nornir_static riff_static mammut_static ${MLPACK_LIBRARY} ${GSL_LIBRARY} pthread) # TODO: Check if dependencies are installed
install(TARGETS nornir_static ARCHIVE DESTINATION lib)

####################
# Uninstall target #
####################
if(NOT TARGET uninstall)
    configure_file(
        "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
        IMMEDIATE @ONLY)

    add_custom_target(uninstall
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake)
endif()
