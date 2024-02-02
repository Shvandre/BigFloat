cmake_minimum_required(VERSION 3.27)
project(bigInt)


set(CMAKE_CXX_STANDARD 20)
Include(FetchContent)

FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

set(CMAKE_CXX_FLAGS "-Ofast")
add_executable(bigInt main.cpp main.h)
target_link_libraries(bigInt PRIVATE Catch2::Catch2WithMain)
