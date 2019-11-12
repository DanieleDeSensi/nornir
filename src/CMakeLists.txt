########################
# Src and header files #
########################
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/nornir/external/fastflow)

file(GLOB SOURCES "*.cpp" "*.c")

if(ENABLE_DATAFLOW)
    file(GLOB SOURCES_DF "dataflow/*.c*")
    list(APPEND SOURCES ${SOURCES_DF})
endif(ENABLE_DATAFLOW)


######################
# Build omp launcher #
######################
if(ENABLE_OMP)
    add_custom_target(omp_launcher
                      COMMAND ${PROJECT_SOURCE_DIR}/bin/generate_omp_launcher.sh ${PROJECT_BINARY_DIR} ${CMAKE_INSTALL_PREFIX})
endif(ENABLE_OMP)

####################
# Library building #
####################

# Must be set before add_library
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

# Download dependencies
include(ExternalProject)

if(ENABLE_OMP)
    set(EXTERNAL_INSTALL_LOCATION_OMP ${PROJECT_BINARY_DIR}/external-omp)
    include_directories(${EXTERNAL_INSTALL_LOCATION_OMP}/include)
    ExternalProject_Add(omp_repo
        GIT_REPOSITORY https://github.com/danieledesensi/LLVM-openmp
        CMAKE_ARGS -DLIBOMP_OMPT_SUPPORT=ON -DLIBOMP_OMPT_OPTIONAL=ON -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION_OMP} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        UPDATE_COMMAND ""
    )
else()
    list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/src/ompt.cpp")
endif(ENABLE_OMP)

if(ENABLE_BLACKBOX)
    set(MAMMUT_PAPI "ON")
else()
    set(MAMMUT_PAPI "OFF")
endif(ENABLE_BLACKBOX)

ExternalProject_Add(mammut_repo
    GIT_REPOSITORY https://github.com/danieledesensi/mammut
    GIT_TAG d290b25
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DENABLE_RAPLCAP=${ENABLE_RAPLCAP} -DENABLE_PAPI=${MAMMUT_PAPI}
    UPDATE_COMMAND ""
)

ExternalProject_Add(riff_repo
    GIT_REPOSITORY https://github.com/danieledesensi/riff
    GIT_TAG 3e1c4f2
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    UPDATE_COMMAND ""
)

if(ENABLE_DATAFLOW)
    # Build external libraries
    add_custom_target(
        compile_ext
        COMMAND make -C ${PROJECT_SOURCE_DIR}/src/external
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
else()
    add_custom_target(
        compile_ext
        COMMAND echo
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif(ENABLE_DATAFLOW)

# Libraries
add_library(nornir SHARED ${SOURCES})
add_library(nornir_static STATIC ${SOURCES})

# Check if mlpack, gsl and armadillo are present
if(ENABLE_GSL)
    find_library(GSL_LIBRARY NAMES gsl)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_GSL")
endif(ENABLE_GSL)

if(ENABLE_MLPACK)
    find_package(MLPACK 3.0.0 REQUIRED)
    include_directories(${MLPACK_INCLUDE_DIRS})
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_MLPACK")

    # Check armadillo deprecated functions
    check_cxx_source_compiles(
    "
    #include <armadillo>
    #include <iostream>
    int main(int argc, char** argv){
        arma::set_cerr_stream(std::cerr);
        return 0;
    }
    " ARMA_DEPRECATED_CERR)

    if(ARMA_DEPRECATED_CERR)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DARMA_DEPRECATED_CERR")
    endif(ARMA_DEPRECATED_CERR)
endif(ENABLE_MLPACK)

if(ENABLE_ARMADILLO)
    find_library(ARMA_LIBRARY NAMES armadillo)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLE_ARMADILLO")
endif(ENABLE_ARMADILLO)

# Find pthreads
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

# TODO: Check if dependencies are installed
if(ENABLE_OMP)
    include(FindOpenMP)
    find_package(OpenMP REQUIRED)
    add_compile_options(OpenMP_CXX_FLAGS)
    set(OMP_DEP "omp_repo")
    set(OMP_DEP_LAUNCHER "omp_launcher")
endif(ENABLE_OMP)

# Dataflow objects
if(ENABLE_DATAFLOW)
    set(DATAFLOW_LIBRARIES ${PROJECT_SOURCE_DIR}/src/external/queues/hzdptr.o ${PROJECT_SOURCE_DIR}/src/external/queues/xxhash.o)
endif(ENABLE_DATAFLOW)

# Dynamic Library
add_dependencies(nornir compile_ext mammut_repo riff_repo ${OMP_DEP} ${OMP_DEP_LAUNCHER})
target_link_libraries(nornir riff mammut ${MLPACK_LIBRARIES} ${GSL_LIBRARY} ${ARMA_LIBRARY} Threads::Threads ${DATAFLOW_LIBRARIES}) 


# Static Library
add_dependencies(nornir_static compile_ext mammut_repo riff_repo ${OMP_DEP} ${OMP_DEP_LAUNCHER})
target_link_libraries(nornir_static riff_static mammut_static  ${MLPACK_LIBRARIES} ${GSL_LIBRARY} ${ARMA_LIBRARY} Threads::Threads ${DATAFLOW_LIBRARIES}) 


###########
# Install #
###########

# Libraries
install(TARGETS nornir LIBRARY DESTINATION lib)
install(TARGETS nornir_static ARCHIVE DESTINATION lib)

# Includes
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/nornir DESTINATION include)

# Install dependencies
install(DIRECTORY ${EXTERNAL_INSTALL_LOCATION}/lib DESTINATION .)
install(DIRECTORY ${EXTERNAL_INSTALL_LOCATION}/include DESTINATION .)

if(ENABLE_OMP)
    install(FILES ${EXTERNAL_INSTALL_LOCATION_OMP}/lib/libomp.so DESTINATION lib RENAME libomp_nornir.so)
endif(ENABLE_OMP)



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
