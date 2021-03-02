find_package(Eigen3 3.3 REQUIRED)

add_library(Reference)
set_property(TARGET Reference PROPERTY CXX_STANDARD 17)

# Globbing may or may not be bad. Use with caution.
file(GLOB_RECURSE Headers include/*.hpp)
target_sources(Reference
    PUBLIC ${Headers})

# Globbing may or may not be bad. Use with caution.
file(GLOB_RECURSE Sources src/*.cpp)
target_sources(Reference
    PRIVATE ${Sources})

target_link_libraries(Reference
    PRIVATE Eigen3::Eigen)

target_include_directories(Reference
    PUBLIC include)
