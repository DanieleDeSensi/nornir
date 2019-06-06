include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_executable(check check.cpp)
target_link_libraries(check LINK_PUBLIC nornir)

add_executable(idlePower idlePower.cpp)
target_link_libraries(idlePower LINK_PUBLIC nornir)

add_executable(ticksPerNs ticksPerNs.cpp)
target_link_libraries(ticksPerNs LINK_PUBLIC nornir)

add_executable(voltageTable voltageTable.cpp)
target_link_libraries(voltageTable LINK_PUBLIC nornir)

add_custom_target(microbench
                  DEPENDS check idlePower ticksPerNs voltageTable
                  COMMAND ${PROJECT_SOURCE_DIR}/microbench/runmicrobenchs_pre.sh ${PROJECT_SOURCE_DIR}
                  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/microbench
                  COMMENT "Analyzing hardware...")
