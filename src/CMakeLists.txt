########################
# Src and header files #
########################
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/nornir/external/fastflow)

file(GLOB SOURCES "*.cpp" "*.c" "dataflow/*.c*")


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

ExternalProject_Add(mammut_repo
    GIT_REPOSITORY https://github.com/danieledesensi/mammut
    GIT_TAG fa14a4b
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    UPDATE_COMMAND ""
)

ExternalProject_Add(riff_repo
    GIT_REPOSITORY https://github.com/danieledesensi/riff
    GIT_TAG 5b341a8
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    UPDATE_COMMAND ""
)

# Build external libraries
add_custom_target(
    compile_ext
    COMMAND make -C ${PROJECT_SOURCE_DIR}/src/external
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
)

# Libraries
add_library(nornir SHARED ${SOURCES})
add_library(nornir_static STATIC ${SOURCES})


# Check if mlpack and gsl are present
find_library(MLPACK_LIBRARY NAMES mlpack)

find_library(GSL_LIBRARY NAMES gsl)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)

# TODO: Check if dependencies are installed
if(ENABLE_OMP)
#    set(OMP_LIBRARY "omp") # Not needed, need to LD_PRELOAD it anyway
    set(OMP_DEP "omp_repo")
    set(OMP_DEP_LAUNCHER "omp_launcher")
endif(ENABLE_OMP)

# Dynamic Library
add_dependencies(nornir compile_ext mammut_repo riff_repo ${OMP_DEP} ${OMP_DEP_LAUNCHER})
target_link_libraries(nornir riff mammut ${MLPACK_LIBRARY} ${GSL_LIBRARY} Threads::Threads ${PROJECT_SOURCE_DIR}/src/external/queues/hzdptr.o ${PROJECT_SOURCE_DIR}/src/external/queues/xxhash.o ${OMP_LIBRARY}) 

# Static Library
add_dependencies(nornir_static compile_ext mammut_repo riff_repo ${OMP_DEP} ${OMP_DEP_LAUNCHER})
target_link_libraries(nornir_static riff_static mammut_static  ${MLPACK_LIBRARY} ${GSL_LIBRARY} Threads::Threads ${PROJECT_SOURCE_DIR}/src/external/queues/hzdptr.o ${PROJECT_SOURCE_DIR}/src/external/queues/xxhash.o ${OMP_LIBRARY}) 


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
    if(EXISTS ${EXTERNAL_INSTALL_LOCATION_OMP}/lib/libomp.so)
        file(RENAME ${EXTERNAL_INSTALL_LOCATION_OMP}/lib/libomp.so ${EXTERNAL_INSTALL_LOCATION_OMP}/lib/libomp_nornir.so)
    endif()
    install(FILES ${EXTERNAL_INSTALL_LOCATION_OMP}/lib/libomp_nornir.so DESTINATION lib)
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
