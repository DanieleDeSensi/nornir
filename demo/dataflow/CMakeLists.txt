include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_executable(asciiMandelbrot asciiMandelbrot.cpp)
target_link_libraries(asciiMandelbrot LINK_PUBLIC nornir)

add_executable(graph graph.cpp)
target_link_libraries(graph LINK_PUBLIC nornir)

add_executable(mapSimple mapSimple.cpp)
target_link_libraries(mapSimple LINK_PUBLIC nornir)

add_executable(pipeline pipeline.cpp)
target_link_libraries(pipeline LINK_PUBLIC nornir)

add_executable(pipeMap2 pipeMap2.cpp)
target_link_libraries(pipeMap2 LINK_PUBLIC nornir)

add_executable(pipeMap pipeMap.cpp)
target_link_libraries(pipeMap LINK_PUBLIC nornir)