cmake_minimum_required(VERSION 3.12.0)
project(mct VERSION 0.2.0 LANGUAGES CXX)

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

option(MCT_ENABLE_MULTITHREAD "Enable multitheading with a posix thread pool" OFF)
if(MCT_ENABLE_MULTITHREAD)
    message(STATUS)
    message(STATUS "Enable multithreading for MCT")
endif()

if(NOT DEFINED MEQ_ENABLE_MULTITHREAD)
    set(MEQ_ENABLE_MULTITHREAD "${MCT_ENABLE_MULTITHREAD}" CACHE BOOL "Forwarded from MCT_ENABLE_MULTITHREAD to MEQ_ENABLE_MULTITHREAD" FORCE)
endif()

add_subdirectory(libmeq)

add_executable(mct src/main.cpp src/Spdata.cpp)
target_include_directories(mct PRIVATE ${PROJECT_SOURCE_DIR}/libmeq)
target_link_libraries(mct PRIVATE meq)

target_compile_definitions(mct PRIVATE $<$<BOOL:${MCT_ENABLE_MULTITHREAD}>:MCT_MULTITHREAD>)

# target_compile_definitions(meq PRIVATE MEQ_ZERNIKE_SERIES_ MEQ_MAX_ZERNIKE_POLAR_ORDER=20)
# target_compile_definitions(mct PRIVATE MEQ_ZERNIKE_SERIES_)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
