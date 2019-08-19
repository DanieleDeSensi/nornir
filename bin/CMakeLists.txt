include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_executable(nornir_blackbox nornir_blackbox.cpp)
#target_include_directories(nornir_blackbox PUBLIC ${PROJECT_SOURCE_DIR}/src/external/tclap-1.2.1/include/)
target_link_libraries(nornir_blackbox LINK_PUBLIC nornir)
install(TARGETS nornir_blackbox DESTINATION bin)

add_executable(nornir_instrumented nornir_instrumented.cpp)
target_link_libraries(nornir_instrumented LINK_PUBLIC nornir)
install(TARGETS nornir_instrumented DESTINATION bin)

add_executable(nornir_instrumented_single nornir_instrumented_single.cpp)
target_link_libraries(nornir_instrumented_single LINK_PUBLIC nornir)
install(TARGETS nornir_instrumented_single DESTINATION bin)

add_executable(nornir_manual_control nornir_manual_control.cpp)
target_include_directories(nornir_manual_control PUBLIC ${PROJECT_SOURCE_DIR}/src/external/tclap-1.2.1/include/)
target_link_libraries(nornir_manual_control LINK_PUBLIC nornir)
install(TARGETS nornir_manual_control DESTINATION bin)

add_subdirectory(nornir_manual_control_web)

if(ENABLE_OMP)
	install(PROGRAMS ${PROJECT_BINARY_DIR}/bin/nornir_openmp DESTINATION bin)
endif(ENABLE_OMP)