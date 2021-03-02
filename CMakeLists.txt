cmake_minimum_required(VERSION 3.17)
project(rewrite LANGUAGES CXX)

add_subdirectory(reference)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Reference)