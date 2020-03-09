cmake_minimum_required(VERSION 3.8..3.14)

if (${CMAKE_VERSION} VERSION_LESS 3.9)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.9)
endif()

# Declare project
project(tqs VERSION 1.0
            DESCRIPTION "Thermal Quench Simulator"
            LANGUAGES CXX)

# Options
option(GIT_SUBMODULE "Check submodules during build" ON)
option(BUILD_TESTS "Build the test framework" OFF)
option(COLOR_TERMINAL "Allow colourful output" ON)

# Add CMake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Get git commit hash
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Generate 'config.h'
configure_file(
    "${PROJECT_SOURCE_DIR}/include/config.h.in"
    "${PROJECT_BINARY_DIR}/include/config.h"
)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory")
endif()

# Fetch SOFTLib
find_package(Git QUIET)
if (GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    if (GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
            WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
            RESULT_VARIABLE GIT_SUBMOD_RESULT)

        if (NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

# Build softlib
add_subdirectory(extern/softlib)

# Sub-directories
add_subdirectory(fvm)
add_subdirectory(src)
add_subdirectory(iface)

#if (BUILD_TESTS)
#    add_subdirectory(tests/c++)
#endif ()

