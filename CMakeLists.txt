cmake_minimum_required(VERSION 3.12.0)
project(mct VERSION 0.1.0)
include(CTest)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebugInfo)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /EHsc")
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra")
endif()

include_directories(.)
add_compile_definitions($<$<CONFIG:Debug>:_DEBUG>)

# add_executable(mct main.cpp)

list(APPEND src_files Contour.cpp GFileRawData.cpp)
list(TRANSFORM src_files PREPEND src/)
list(APPEND tests vec-test util-test read-test)

list(LENGTH tests test_num)
message(STATUS)
message(STATUS  ${test_num} " tests found: ")
enable_testing()
foreach(test ${tests})
    string(JOIN "" test_src_file test/ ${test} .cpp)
    add_executable(${test} ${test_src_file} ${src_files})

    string(JOIN "-" test_name mct ${test})
    add_test(${test_name} ${test})
    message(STATUS ${test} " from source file " ${test_src_file})
endforeach(test ${tests})
message(STATUS)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
