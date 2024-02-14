cmake_minimum_required(VERSION 3.11)
project(BigFloat
        VERSION 1.0.0
        DESCRIPTION "Library that implements fixed point numbers"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "-O0 -g")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS "-Ofast -pthread")
else()
    message("Unknown build type")
endif()

#Just build the library target
add_library(BigFloat BigFloat.cpp BigFloat.h)

#Test target
Include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.4.0 # or a later release
)
FetchContent_MakeAvailable(Catch2)

add_executable(tests tests.cpp BigFloat.cpp BigFloat.h)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain)

#Calc Pi tagret
add_executable(CalcPi main.cpp BigFloat)