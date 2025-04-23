cmake_minimum_required(VERSION 3.10)
# set(CMAKE_C_COMPILER "/usr/bin/clang")
# set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
set(CMAKE_TRY_COMPILE_TARGET_TYPE "STATIC_LIBRARY")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
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

# vector - load
add_executable(astar-maskedgather "astar.cpp" "util.h")
target_compile_options(astar-maskedgather PUBLIC "-mavx2")
target_compile_definitions(astar-maskedgather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDGATHER)
target_link_libraries(astar-maskedgather pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown64 "astar.cpp" "util.h")
target_compile_options(astar-rounddown64 PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown64 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_64)
target_link_libraries(astar-rounddown64 pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown32 "astar.cpp" "util.h")
target_compile_options(astar-rounddown32 PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown32 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_32)
target_link_libraries(astar-rounddown32 pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown16 "astar.cpp" "util.h")
target_compile_options(astar-rounddown16 PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown16 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_16)
target_link_libraries(astar-rounddown16 pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown8 "astar.cpp" "util.h")
target_compile_options(astar-rounddown8 PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown8 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_8)
target_link_libraries(astar-rounddown8 pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast64 "astar.cpp" "util.h")
target_compile_options(astar-maskedlast64 PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast64 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_64)
target_link_libraries(astar-maskedlast64 pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast32 "astar.cpp" "util.h")
target_compile_options(astar-maskedlast32 PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast32 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_32)
target_link_libraries(astar-maskedlast32 pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast16 "astar.cpp" "util.h")
target_compile_options(astar-maskedlast16 PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast16 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_16)
target_link_libraries(astar-maskedlast16 pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast8 "astar.cpp" "util.h")
target_compile_options(astar-maskedlast8 PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast8 PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_8)
target_link_libraries(astar-maskedlast8 pthread ${Boost_LIBRARIES})


# vector - gather
add_executable(astar-maskedgather-gather "astar.cpp" "util.h")
target_compile_options(astar-maskedgather-gather PUBLIC "-mavx2")
target_compile_definitions(astar-maskedgather-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDGATHER AWU_VEC_GATHER)
target_link_libraries(astar-maskedgather-gather pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown64-gather "astar.cpp" "util.h")
target_compile_options(astar-rounddown64-gather PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown64-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_64 AWU_VEC_GATHER)
target_link_libraries(astar-rounddown64-gather pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown32-gather "astar.cpp" "util.h")
target_compile_options(astar-rounddown32-gather PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown32-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_32 AWU_VEC_GATHER)
target_link_libraries(astar-rounddown32-gather pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown16-gather "astar.cpp" "util.h")
target_compile_options(astar-rounddown16-gather PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown16-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_16 AWU_VEC_GATHER)
target_link_libraries(astar-rounddown16-gather pthread ${Boost_LIBRARIES})

add_executable(astar-rounddown8-gather "astar.cpp" "util.h")
target_compile_options(astar-rounddown8-gather PUBLIC "-mavx2")
target_compile_definitions(astar-rounddown8-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_ROUNDDOWN AWU_MASKWIDTH_8 AWU_VEC_GATHER)
target_link_libraries(astar-rounddown8-gather pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast64-gather "astar.cpp" "util.h")
target_compile_options(astar-maskedlast64-gather PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast64-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_64 AWU_VEC_GATHER)
target_link_libraries(astar-maskedlast64-gather pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast32-gather "astar.cpp" "util.h")
target_compile_options(astar-maskedlast32-gather PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast32-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_32 AWU_VEC_GATHER)
target_link_libraries(astar-maskedlast32-gather pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast16-gather "astar.cpp" "util.h")
target_compile_options(astar-maskedlast16-gather PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast16-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_16 AWU_VEC_GATHER)
target_link_libraries(astar-maskedlast16-gather pthread ${Boost_LIBRARIES})

add_executable(astar-maskedlast8-gather "astar.cpp" "util.h")
target_compile_options(astar-maskedlast8-gather PUBLIC "-mavx2")
target_compile_definitions(astar-maskedlast8-gather PUBLIC AWU_MASK_VECTOR AWU_VECSCHEME_MASKEDLAST AWU_MASKWIDTH_8 AWU_VEC_GATHER)
target_link_libraries(astar-maskedlast8-gather pthread ${Boost_LIBRARIES})

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