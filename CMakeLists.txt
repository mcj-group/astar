cmake_minimum_required(VERSION 3.10)
project(astar)

include_directories(${Boost_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../CPS/include)

add_executable(astar "astar.cpp" "util.h")
target_compile_options(astar PUBLIC "-mavx2")
target_link_libraries(astar pthread ${Boost_LIBRARIES})
