add_executable(example1 example1.cpp)
target_link_libraries(example1 PRIVATE reference)

add_custom_command(TARGET example1 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_LIST_DIR}/data $<TARGET_FILE_DIR:example1>/data)