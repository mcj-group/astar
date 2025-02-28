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

add_executable(astar-scal "astar.cpp" "util.h")
target_compile_options(astar-scal PUBLIC "-mavx2")
target_compile_definitions(astar-scal PUBLIC AWU_MASK_SCALAR)
target_link_libraries(astar-scal pthread ${Boost_LIBRARIES})

add_executable(astar-valways "astar.cpp" "util.h")
target_compile_options(astar-valways PUBLIC "-mavx2")
target_compile_definitions(astar-valways PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ALWAYS)
target_link_libraries(astar-valways pthread ${Boost_LIBRARIES})

add_executable(astar-vvecfill32 "astar.cpp" "util.h")
target_compile_options(astar-vvecfill32 PUBLIC "-mavx2")
target_compile_definitions(astar-vvecfill32 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_FILLVEC)
target_link_libraries(astar-vvecfill32 pthread ${Boost_LIBRARIES})

add_executable(astar-vvecfill64 "astar.cpp" "util.h")
target_compile_options(astar-vvecfill64 PUBLIC "-mavx2")
target_compile_definitions(astar-vvecfill64 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_FILLVEC AWU_MASKWIDTH_64)
target_link_libraries(astar-vvecfill64 pthread ${Boost_LIBRARIES})

add_executable(astar-vmov32 "astar.cpp" "util.h")
target_compile_options(astar-vmov32 PUBLIC "-mavx2")
target_compile_definitions(astar-vmov32 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MOSTOFVEC)
target_link_libraries(astar-vmov32 pthread ${Boost_LIBRARIES})

add_executable(astar-vmov64 "astar.cpp" "util.h")
target_compile_options(astar-vmov64 PUBLIC "-mavx2")
target_compile_definitions(astar-vmov64 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MOSTOFVEC AWU_MASKWIDTH_64)
target_link_libraries(astar-vmov64 pthread ${Boost_LIBRARIES})

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