include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_executable(manager-blackbox manager-blackbox.cpp)
#target_include_directories(manager-blackbox PUBLIC ${PROJECT_SOURCE_DIR}/src/external/tclap-1.2.1/include/)
target_link_libraries(manager-blackbox LINK_PUBLIC nornir)

add_executable(manager-external manager-external.cpp)
target_link_libraries(manager-external LINK_PUBLIC nornir)

add_executable(manager-single manager-single.cpp)
target_link_libraries(manager-single LINK_PUBLIC nornir)

add_executable(selector-manual selector-manual.cpp)
target_include_directories(selector-manual PUBLIC ${PROJECT_SOURCE_DIR}/src/external/tclap-1.2.1/include/)
target_link_libraries(selector-manual LINK_PUBLIC nornir)

add_subdirectory(selector-manual-web)