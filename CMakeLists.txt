cmake_minimum_required(VERSION 3.24)

project(c2c-sepia)

if (UNIX)
  set (CMAKE_CXX_FLAGS "-g3 -Wall -O0")
endif ()

set (CMAKE_BINARY_DIR "build")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

add_executable(c2c-sepia src/main.cxx)