cmake_minimum_required(VERSION 3.17)
project(rewrite LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(reference)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Reference)

set(Resources 1ubq.pdb)
foreach(resource IN LISTS Resources)
    configure_file(${resource} ${resource} COPYONLY)
endforeach()
