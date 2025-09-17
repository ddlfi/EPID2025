cmake_minimum_required(VERSION 3.10)
project(SignatureProject)

set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)

# 设置构建类型为Release（如果未指定）
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# 激进的优化标志
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -flto")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native -flto")

# 链接时优化
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "-flto")

find_package(OpenSSL REQUIRED)
find_package(OpenMP REQUIRED)

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-fcommon)
endif()

# 向量化和加密指令优化
add_compile_options(-mfma -mpclmul -msse2 -mavx2 -maes)

add_subdirectory(sha3)

include_directories(${CMAKE_SOURCE_DIR}/include)

file(GLOB_RECURSE PROJECT_SOURCES "src/*.c" "src/*.cpp")

add_executable(cases_test ${PROJECT_SOURCES})

target_link_libraries(cases_test PRIVATE sha3)
target_link_libraries(cases_test PRIVATE OpenSSL::Crypto)
target_link_libraries(cases_test PRIVATE OpenMP::OpenMP_CXX)
