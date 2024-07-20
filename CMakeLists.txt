cmake_minimum_required(VERSION 3.1)
project(ChemPlasKin)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(Threads REQUIRED)
find_library(ACCELERATE_FRAMEWORK Accelerate)

include_directories("../cantera/include")
include_directories("../cantera/src")
include_directories("../cantera/include/cantera/ext")
link_directories("../cantera/build/lib")

add_executable(ChemPlasKin src/master.cpp)

target_link_libraries(ChemPlasKin cantera_shared ${ACCELERATE_FRAMEWORK} Threads::Threads)

# Set the RPATH to use the relative path
set_target_properties(ChemPlasKin PROPERTIES
    BUILD_WITH_INSTALL_RPATH TRUE
    INSTALL_RPATH "@loader_path/../../cantera/build/lib"
    SKIP_BUILD_RPATH FALSE
    BUILD_RPATH "../../cantera/build/lib"
)


