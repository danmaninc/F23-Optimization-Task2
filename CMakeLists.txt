cmake_minimum_required(VERSION 3.23)
project(optimization)

set(CMAKE_CXX_STANDARD 23)

add_executable(optimization main.cpp)
target_link_libraries(optimization gtest gmock pthread)