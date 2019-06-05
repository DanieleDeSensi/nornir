include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

add_executable(application_instrumented application_instrumented.cpp)
target_link_libraries(application_instrumented LINK_PUBLIC nornir)

add_executable(farm_accelerator_1 farm_accelerator_1.cpp)
target_link_libraries(farm_accelerator_1 LINK_PUBLIC nornir)

add_executable(farm_accelerator_2 farm_accelerator_2.cpp)
target_link_libraries(farm_accelerator_2 LINK_PUBLIC nornir)

add_executable(farm_accelerator_3 farm_accelerator_3.cpp)
target_link_libraries(farm_accelerator_3 LINK_PUBLIC nornir)

add_executable(farm_accelerator_4 farm_accelerator_4.cpp)
target_link_libraries(farm_accelerator_4 LINK_PUBLIC nornir)

add_executable(farm farm.cpp)
target_link_libraries(farm LINK_PUBLIC nornir)

add_executable(farm_ffcompat farm_ffcompat.cpp)
target_link_libraries(farm_ffcompat LINK_PUBLIC nornir)
target_include_directories(farm_ffcompat PUBLIC ${PROJECT_SOURCE_DIR}/include/nornir/external/fastflow/)

add_executable(farm_ffcompat_new_api farm_ffcompat_new_api.cpp)
target_link_libraries(farm_ffcompat_new_api LINK_PUBLIC nornir)
target_include_directories(farm_ffcompat_new_api PUBLIC ${PROJECT_SOURCE_DIR}/include/nornir/external/fastflow/)

add_executable(parallel_for parallel_for.cpp)
target_link_libraries(parallel_for LINK_PUBLIC nornir)
target_include_directories(parallel_for PUBLIC ${PROJECT_SOURCE_DIR}/include/nornir/external/fastflow/)

add_executable(pipe_ffcompat pipe_ffcompat.cpp)
target_link_libraries(pipe_ffcompat LINK_PUBLIC nornir)
target_include_directories(pipe_ffcompat PUBLIC ${PROJECT_SOURCE_DIR}/include/nornir/external/fastflow/)