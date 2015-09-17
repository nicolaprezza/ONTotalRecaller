cmake_minimum_required(VERSION 2.6)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

#set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/.. )

project (ONTotalRecaller)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/include/internal)
include_directories(${PROJECT_SOURCE_DIR}/extern/fast5/src)
include_directories(/usr/local/hdf5/include)
include_directories(~/include)

LINK_DIRECTORIES(/usr/local/hdf5/lib)
LINK_DIRECTORIES(~/lib)

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "--std=c++11")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g")
set(CMAKE_CXX_FLAGS_RELEASE "-g -ggdb -Ofast -fstrict-aliasing -DNDEBUG -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

add_executable(ontrc_param ontrc_param.cpp)
TARGET_LINK_LIBRARIES(ontrc_param hdf5)
TARGET_LINK_LIBRARIES(ontrc_param sdsl)
TARGET_LINK_LIBRARIES(ontrc_param divsufsort)
TARGET_LINK_LIBRARIES(ontrc_param divsufsort64)

add_executable(ontrc_build ontrc_build.cpp)
TARGET_LINK_LIBRARIES(ontrc_build hdf5)
TARGET_LINK_LIBRARIES(ontrc_build sdsl)
TARGET_LINK_LIBRARIES(ontrc_build divsufsort)
TARGET_LINK_LIBRARIES(ontrc_build divsufsort64)

add_executable(ontrc_call ontrc_call.cpp)
TARGET_LINK_LIBRARIES(ontrc_call hdf5)
TARGET_LINK_LIBRARIES(ontrc_call sdsl)
TARGET_LINK_LIBRARIES(ontrc_call divsufsort)
TARGET_LINK_LIBRARIES(ontrc_call divsufsort64)

add_executable(debug debug.cpp)
TARGET_LINK_LIBRARIES(debug hdf5)
TARGET_LINK_LIBRARIES(debug sdsl)
TARGET_LINK_LIBRARIES(debug divsufsort)
TARGET_LINK_LIBRARIES(debug divsufsort64)

add_executable(ontrc_test ontrc_test.cpp)
TARGET_LINK_LIBRARIES(ontrc_test hdf5)
TARGET_LINK_LIBRARIES(ontrc_test sdsl)
TARGET_LINK_LIBRARIES(ontrc_test divsufsort)
TARGET_LINK_LIBRARIES(ontrc_test divsufsort64)