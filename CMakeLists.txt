cmake_minimum_required(VERSION 2.8)
project(rns-heaan CXX ASM-ATT)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include(CMakeCompilers.txt)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
include_directories(${GMP_INCLUDES} ${MPFR_INCLUDES} /usr/local/include)
link_directories(/usr/local/lib)
set(RUN_LIB yell boost_system ${MPFR_LIBRARIES} ${GMP_LIBRARIES})

# add_definitions(-DYELL_USE_AVX)
#add_definitions(-DYELL_USE_32BITS_MODULI)
#add_definitions(-DYELL_USE_MEM_POOL)
add_subdirectory(yell)
include_directories(include/ yell/include)

file(GLOB fHE_SRC src/*.cpp)
add_library(fHE STATIC ${fHE_SRC})

add_executable(main main.cpp)
target_link_libraries(main fHE ${RUN_LIB})

find_package(GTest)
# if(GTEST_FOUND)
#   enable_testing()
#   add_subdirectory(test)
# endif(GTEST_FOUND)
