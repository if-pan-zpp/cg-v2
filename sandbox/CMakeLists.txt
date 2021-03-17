add_executable(main main.cpp)
target_link_libraries(main TARGET reference)

add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_LIST_DIR}/data $<TARGET_FILE_DIR:main>)