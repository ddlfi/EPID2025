cmake_minimum_required(VERSION 3.10)
project(SignatureProject)


set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
find_package(OpenSSL REQUIRED)

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-fcommon)
endif()

set(CMAKE_CXX_FLAGS "-std=c++17 -mfma")
add_compile_options(-mfma)

add_subdirectory(sha3)


include_directories(${CMAKE_SOURCE_DIR}/include)


file(GLOB_RECURSE PROJECT_SOURCES "src/*.c" "src/*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpclmul -msse2")


add_executable(cases_test ${PROJECT_SOURCES})


target_compile_options(cases_test PRIVATE -mavx2)


target_link_libraries(cases_test PRIVATE sha3)
target_link_libraries(cases_test PRIVATE OpenSSL::Crypto)