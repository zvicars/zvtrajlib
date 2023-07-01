cmake_minimum_required(VERSION 3.10)

project(sphericart LANGUAGES C CXX)

file(READ ${PROJECT_SOURCE_DIR}/VERSION SPHERICART_VERSION)
string(STRIP ${SPHERICART_VERSION} SPHERICART_VERSION)
string(REGEX REPLACE "^([0-9]+)\\..*" "\\1" SPHERICART_VERSION_MAJOR "${SPHERICART_VERSION}")
string(REGEX REPLACE "^[0-9]+\\.([0-9]+).*" "\\1" SPHERICART_VERSION_MINOR "${SPHERICART_VERSION}")
string(REGEX REPLACE "^[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" SPHERICART_VERSION_PATCH "${SPHERICART_VERSION}")

OPTION(SPHERICART_BUILD_TESTS "Build and run tests for Sphericart" OFF)
OPTION(SPHERICART_OPENMP "Try to use OpenMP when compiling Sphericart" ON)
OPTION(SPHERICART_ARCH_NATIVE "Try to use -march=native when compiling Sphericart" ON)


add_library(sphericart SHARED
    "src/sphericart.cpp"
    "include/sphericart.hpp"
)

set_target_properties(sphericart PROPERTIES
    VERSION ${SPHERICART_VERSION}
    SOVERSION ${SPHERICART_VERSION_MAJOR}.${SPHERICART_VERSION_MINOR}
    POSITION_INDEPENDENT_CODE ON
)

# we need to compile sphericart with C++17 for if constexpr
target_compile_features(sphericart PRIVATE cxx_std_17)

target_include_directories(sphericart PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include )
# Handle optimization and OpenMP flags
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wunknown-pragmas" COMPILER_SUPPORTS_WPRAGMAS)
if (SPHERICART_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP is enabled")
        target_link_libraries(sphericart PUBLIC OpenMP::OpenMP_CXX)
    else()
        message(WARNING "Could not find OpenMP")
        if(COMPILER_SUPPORTS_WPRAGMAS)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
        endif()
    endif()
else()
    if(COMPILER_SUPPORTS_WPRAGMAS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    endif()
endif()