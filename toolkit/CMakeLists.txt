find_package(Eigen3 3.3 REQUIRED)

set(Target base)
add_library(${Target})

# Globbing may or may not be bad. Use with caution.
file(GLOB_RECURSE Headers include/*.hpp)
target_sources(${Target}
    PUBLIC ${Headers})

# Globbing may or may not be bad. Use with caution.
file(GLOB_RECURSE Sources src/*.cpp)
target_sources(${Target}
    PRIVATE ${Sources})

target_link_libraries(${Target}
    PUBLIC Eigen3::Eigen)

target_include_directories(${Target}
    PUBLIC include)
