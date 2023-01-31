cmake_minimum_required(VERSION 3.24)

project(c2c-sepia)

if (UNIX)
  #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -lblas -llapack")   #blas is very fast for me on this machine
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -lopenblas") 
endif ()

add_definitions(-DARMA_DONT_USE_WRAPPER)
set (CMAKE_BINARY_DIR "build")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(EXEC_MAIN "src/main.cxx")
set(EXPERIMENT_MAIN "src/testingPerformance.cxx")
set(LIB_SRCS 
    src/WeightedEdgeGraph.cxx
    src/utilities.cxx
    src/optimization.cxx
    src/Matrix.cxx
    src/Computation.cxx
)

if(CUDAavailable)
#using CUDA implementation
elseif(OpenCLavailable)
#using OpenCL implementation
elseif(REPASTHPCavailable)
#using Repast HPC implementation
endif(CUDAavailable)


find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_library( mysharedlib SHARED ${LIB_SRCS} )

add_executable( c2c-sepia ${EXEC_MAIN} )
target_link_libraries( c2c-sepia mysharedlib )

add_executable( experiments ${EXPERIMENT_MAIN})
target_link_libraries( experiments mysharedlib )


# TESTING
enable_testing()

add_executable(GraphTesting  "src/GraphTesting.cc")
target_link_libraries(
  GraphTesting
  GTest::gtest_main
  mysharedlib
)


add_executable(MatrixTesting  "src/MatrixTesting.cc")
target_link_libraries(
  MatrixTesting
  GTest::gtest_main
  mysharedlib
)

add_executable(ComputationTesting  "src/ComputationTesting.cc")
target_link_libraries(
  ComputationTesting
  GTest::gtest_main
  mysharedlib
)

include(GoogleTest)
gtest_discover_tests(GraphTesting)
gtest_discover_tests(MatrixTesting)
gtest_discover_tests(ComputationTesting)