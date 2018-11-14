cmake_minimum_required(VERSION 3.10)
project(code)

set(CMAKE_CXX_STANDARD 11)

add_executable(code testMain.cpp Tree.h Tree.cpp try.cpp main.cpp testRand.cpp)