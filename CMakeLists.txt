cmake_minimum_required (VERSION 3.8)
project (Spatial_context)

set(CMAKE_BUILD_TYPE "Release")
# set(CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -g -pthread")


find_package(PCL REQUIRED)
find_package(OpenCV REQUIRED)
find_package ( OpenMP )
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


include_directories(${PCL_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
include_directories("./Spatialcontext/")

add_executable(spatialcontext main.cpp  Spatialcontext/Spatialcontext.cpp Spatialcontext/nanoflann.hpp)
target_link_libraries(spatialcontext  OpenMP::OpenMP_CXX ${PCL_LIBRARIES} ${OpenCV_LIBS})