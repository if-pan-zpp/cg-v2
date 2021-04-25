set(TARGET mdk)

add_library(${TARGET})

file(GLOB_RECURSE PUBLIC_SOURCES include/*)
file(GLOB_RECURSE PRIVATE_SOURCES src/*)

target_sources(${TARGET}
    PUBLIC
        ${PUBLIC_SOURCES}
    PRIVATE
    ${PRIVATE_SOURCES})

target_include_directories(${TARGET}
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${TARGET}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include/${TARGET}>
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )

target_compile_options(${TARGET}
    PRIVATE
        -Wall -Wextra -fPIC
        $<$<CONFIG:Debug>:-pg>
        $<$<CONFIG:Release>:-Ofast -march=native>
    )

target_compile_features(${TARGET}
    PRIVATE
        cxx_std_17)

#add_subdirectory(vendor)

find_package(Eigen3 3.3 REQUIRED)
target_link_libraries(${TARGET} PUBLIC
    Eigen3::Eigen)
