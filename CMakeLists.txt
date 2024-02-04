cmake_minimum_required(VERSION 3.27)
project(BigFloat
        VERSION 1.0.0
        DESCRIPTION "Library that implements fixed point numbers"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

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


set(CMAKE_CXX_FLAGS "-Ofast")

#Main tagret
add_executable(CalcPi main.cpp BigFloat.cpp BigFloat.h)