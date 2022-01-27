add_executable(test1 main.cpp)
target_link_libraries(test1 PRIVATE reference)

add_custom_command(TARGET test1 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_LIST_DIR}/data $<TARGET_FILE_DIR:test1>/data)
