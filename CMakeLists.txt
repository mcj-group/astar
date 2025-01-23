cmake_minimum_required(VERSION 3.10)
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
project(astar)

include_directories(${Boost_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../CPS/include)

add_executable(astar "astar.cpp" "util.h")
target_compile_options(astar PUBLIC "-mavx2")
target_link_libraries(astar pthread ${Boost_LIBRARIES})

# # Define the gcov target
# add_custom_target(gcov
#     COMMAND ${CMAKE_COMMAND} -E echo "Building with GCOV enabled"
#     COMMAND ${CMAKE_COMMAND} -E make_directory gcov_build
#     COMMAND ${CMAKE_COMMAND}
#         -DCMAKE_BUILD_TYPE=Debug
#         -DCMAKE_C_FLAGS="-g --coverage"
#         -DCMAKE_CXX_FLAGS="-g --coverage"
#         -DCMAKE_EXE_LINKER_FLAGS="--coverage"
#         -Bgcov_build -H.
#     COMMAND ${CMAKE_COMMAND} --build gcov_build --target astar
#     WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
#     COMMENT "Building with coverage flags"
# )