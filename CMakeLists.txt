cmake_minimum_required(VERSION 3.12..3.14)

if (${CMAKE_VERSION} VERSION_LESS 3.13)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
    cmake_policy(VERSION 3.13)
endif()

# Ensure that we build in optimized (Release) mode
# if not explicitly stated otherwise
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(CMakeDependentOption)

# Declare project
project(dream VERSION 1.0
            DESCRIPTION "Disruption and Runaway Electron Avoidance Model"
            LANGUAGES CXX)

# Options
option(COLOR_TERMINAL "Allow colourful output" ON)
option(DREAM_BUILD_TESTS "Build the test framework" ON)
option(DREAM_BUILD_PYFACE "Build the DREAM Python interface" OFF)
option(GIT_SUBMODULE "Check submodules during build" ON)
option(PETSC_WITH_MPI "If ON, indicates that PETSc was linked with MPI" ON)

# Require position-indepdent code when building PyFace
cmake_dependent_option(
    POS_INDP_CODE "Generate position-independent code in static libraries" OFF
    "DREAM_BUILD_PYFACAE" ON
)

# Add CMake modules
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Get git commit hash
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Generate 'config.h'
configure_file(
    "${PROJECT_SOURCE_DIR}/include/config.h.in"
    "${PROJECT_BINARY_DIR}/include/FVM/config.h"
)
configure_file(
    "${PROJECT_SOURCE_DIR}/include/config.h.in"
    "${PROJECT_BINARY_DIR}/include/tests/cxx/config.h"
)
configure_file(
    "${PROJECT_SOURCE_DIR}/include/config.h.in"
    "${PROJECT_BINARY_DIR}/include/DREAM/config.h"
)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory")
endif()

# Find git
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

# Enable compilation warnings on all targets
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-ignored-qualifiers)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic -Wformat-security -Wno-misleading-indentation -Wno-ignored-qualifiers)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-misleading-indentation -Wno-ignored-qualifiers)
    if (CMAKE_BUILD_TYPE STREQUAL "Release")
        add_compile_options(-O3)
    endif()
else ()
    message(WARNING "You are using an un-tested compiler.")
endif ()


# Build softlib
add_subdirectory(extern/softlib)

# Sub-directories
add_subdirectory(fvm)
add_subdirectory(src)
add_subdirectory(iface)

if (DREAM_BUILD_PYFACE)
    add_subdirectory(pyface)
endif()

if (DREAM_BUILD_TESTS)
    add_subdirectory(tests/cxx)
endif ()

