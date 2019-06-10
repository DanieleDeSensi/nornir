# additional target to perform cppcheck run, requires cppcheck
find_package(cppcheck)

if(!CPPCHECK_FOUND)
  message("cppcheck not found. Please install it to run 'make cppcheck'")
endif()

add_custom_target(
        cppcheck
        COMMAND ${CPPCHECK_EXECUTABLE} --xml --xml-version=2 --inline-suppr --enable=warning,performance,information,style --inline-suppr --error-exitcode=1 --force  ${PROJECT_SOURCE_DIR} --suppressions-list=${PROJECT_SOURCE_DIR}/test/cppcheck/suppressions-list.txt -i${PROJECT_SOURCE_DIR}/src/external/queues -i${PROJECT_BINARY_DIR} -i${PROJECT_SOURCE_DIR}/test -i${PROJECT_SOURCE_DIR}/cmake -i${PROJECT_SOURCE_DIR}/demo  2> ${PROJECT_BINARY_DIR}/cppcheck-report.xml || (cat ${PROJECT_BINARY_DIR}/cppcheck-report.xml; exit 2) 
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
)


