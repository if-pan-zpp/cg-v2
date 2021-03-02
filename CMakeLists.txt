cmake_minimum_required(VERSION 3.17)
project(rewrite LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)

add_subdirectory(reference)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE Reference)

add_custom_command(
    TARGET main PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_LIST_DIR}/data $<TARGET_FILE_DIR:main>/data)
