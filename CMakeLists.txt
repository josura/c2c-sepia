cmake_minimum_required(VERSION 3.24)

project(c2c-sepia)

#if (UNIX)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fopenmp")   #  -lblas -llapack are included directly #blas is very fast for me on this machine
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -lopenblas")

add_definitions(-DARMA_DONT_USE_WRAPPER)
set (CMAKE_BINARY_DIR "build")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
else()
  message( SEND_ERROR "lapack or blas not found")
endif()


FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

set(EXEC_MAIN "src/main.cxx")
set(EXPERIMENT_MAIN "src/testingPerformance.cxx")
set(LIB_SRCS 
    src/WeightedEdgeGraph.cxx
    src/utilities.cxx
    src/armaUtilities.cxx
    src/optimization.cxx
    src/Matrix.cxx
    src/Computation.cxx
    src/DissipationModel.cxx
    src/DissipationModelPow.cxx
    src/DissipationModelRandom.cxx
    src/DissipationModelPeriodic.cxx
    src/DissipationModelScaled.cxx
)

if(CUDAavailable)
#using CUDA implementation
elseif(OpenCLavailable)
#using OpenCL implementation
elseif(REPASTHPCavailable)
#using Repast HPC implementation
endif(CUDAavailable)

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
target_link_libraries( c2c-sepia 
                       mysharedlib
                       ${lapackblas_libraries}
                       ${Boost_LIBRARIES})

add_executable( experiments ${EXPERIMENT_MAIN})
target_link_libraries( experiments 
                       mysharedlib  
                       ${lapackblas_libraries})


# TESTING
enable_testing()

add_executable(GraphTesting  "src/GraphTesting.cc")
target_link_libraries(
  GraphTesting
  GTest::gtest_main
  mysharedlib
  ${lapackblas_libraries}
)


add_executable(MatrixTesting  "src/MatrixTesting.cc")
target_link_libraries(
  MatrixTesting
  GTest::gtest_main
  mysharedlib
  ${lapackblas_libraries}
)

add_executable(ComputationTesting  "src/ComputationTesting.cc")
target_link_libraries(
  ComputationTesting
  GTest::gtest_main
  mysharedlib
  ${lapackblas_libraries}
)

add_executable(DissipationModelTesting  "src/DissipationModelTesting.cc")
target_link_libraries(
  DissipationModelTesting
  GTest::gtest_main
  mysharedlib
  ${lapackblas_libraries}
)

include(GoogleTest)
gtest_discover_tests(GraphTesting)
gtest_discover_tests(MatrixTesting)
gtest_discover_tests(ComputationTesting)
gtest_discover_tests(DissipationModelTesting)