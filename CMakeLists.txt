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

add_subdirectory(lib)

add_executable(mct src/main.cpp src/Spdata.cpp)
target_include_directories(mct PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(mct PRIVATE meq)

target_compile_definitions(meq PRIVATE MEQ_ZERNIKE_SERIES_ MEQ_MAX_ZERNIKE_POLAR_ORDER=20)
target_compile_definitions(mct PRIVATE MEQ_ZERNIKE_SERIES_)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
