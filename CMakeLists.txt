cmake_minimum_required(VERSION 3.17)
project(rewrite LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

option(LEGACY_MODE "Work as close to f77 version as possible." OFF)
if (LEGACY_MODE)
  add_definitions(-DLEGACY_MODE)
endif()

add_subdirectory(toolkit)
add_subdirectory(reference)
add_subdirectory(sandbox)
add_subdirectory(tests)
