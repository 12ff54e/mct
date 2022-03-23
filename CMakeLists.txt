cmake_minimum_required(VERSION 3.0.0)
project(mct VERSION 0.1.0)

add_compile_options(-std=c++11)
add_compile_options(-Wall)
include_directories($ENV{HOME}/software/install/include)
include_directories(.)

# add_executable(mct main.cpp)
add_executable(util_test test/util_test.cpp)
add_executable(vec_test test/vec_test.cpp)
add_executable(read_test test/read_test.cpp src/GFileRawData.cpp)

include(CTest)
enable_testing()
add_test(util util_test)
add_test(vec vec_test)
add_test(read-gfile read_test)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
