cmake_minimum_required(VERSION 3.12.0)
project(meq VERSION 1.0.0 LANGUAGES CXX)
include(CTest)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebugInfo)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    set(CMAKE_CXX_FLAGS "/W4 /EHsc")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-pre-c++14-compat -Wno-padded -Wno-exit-time-destructors -Wno-unsafe-buffer-usage -Wno-covered-switch-default")
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra")
endif()

add_compile_definitions($<$<OR:$<CONFIG:Debug>,$<CONFIG:RelWithDebInfo>>:_DEBUG>)

list(APPEND lib_src_files Contour.cpp GFileRawData.cpp MagneticEquilibrium.cpp)
add_library(meq ${lib_src_files})

list(APPEND tests vec-test util-test poly-test bigint-test zernike-test)
list(LENGTH tests test_num)
message(STATUS)
message(STATUS  ${test_num} " tests found: ")
enable_testing()
foreach(test ${tests})
    string(JOIN "" test_src_file test/ ${test} .cpp)
    add_executable(${test} ${test_src_file})
    target_include_directories(${test} PRIVATE ${PROJECT_SOURCE_DIR})

    string(JOIN "-" test_name mct ${test})
    add_test(${test_name} ${test})
    message(STATUS ${test} " from source file " ${test_src_file})
endforeach(test ${tests})
message(STATUS)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
