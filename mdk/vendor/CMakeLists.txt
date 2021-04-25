include(FetchContent)

# This is for suppressing warnings from Eigen3
set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")

message("Fetching: Eigen3")
FetchContent_Declare(
    Eigen3
    URL https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.zip
)
FetchContent_MakeAvailable(Eigen3)

message("Fetching: Catch2")
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v2.13.1)
FetchContent_MakeAvailable(Catch2)

message("Fetching done")