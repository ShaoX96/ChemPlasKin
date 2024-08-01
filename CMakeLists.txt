cmake_minimum_required(VERSION 3.1)
project(ChemPlasKin)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(Threads REQUIRED)

include_directories("../cantera/include")
include_directories("../cantera/src")
include_directories("../cantera/include/cantera/ext")
link_directories("../cantera/build/lib")

add_executable(ChemPlasKin src/master.cpp)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    target_link_libraries(ChemPlasKin PRIVATE cantera_shared ${ACCELERATE_FRAMEWORK} Threads::Threads)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_libraries(ChemPlasKin PRIVATE cantera_shared stdc++fs Threads::Threads)
endif()
