cmake_minimum_required(VERSION 3.24)

project(c2c-sepia)

if (UNIX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 ")
endif ()

set (CMAKE_BINARY_DIR "build")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(EXEC_MAIN "src/main.cxx")
#set(EXPERIMENT_MAIN "experiments.cxx")
set(LIB_SRCS 
    src/WeightedEdgeGraph.cxx
    src/utilities.cxx
)

add_library( mysharedlib SHARED ${LIB_SRCS} )

add_executable( c2c-sepia ${EXEC_MAIN} )
target_link_libraries( c2c-sepia mysharedlib )

#add_executable( experiments ${EXPERIMENT_MAIN})
#target_link_libraries( experiments mysharedlib )
