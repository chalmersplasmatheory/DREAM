
set(tqstests_core
    "${PROJECT_SOURCE_DIR}/tests/cxx/runtest.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/UnitTest.cpp"
)

set(tqstests_fvm
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/RadialGrid.cpp"
)

add_executable(tqstests ${tqstests_core} ${tqstests_fvm})
target_include_directories(tqstests PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/tests/cxx/include")
target_link_libraries(tqstests PUBLIC fvm)
target_link_libraries(tqstests PUBLIC tqs)

# Require C++17
set_target_properties(tqstests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

######################
# MODULES
######################
# Math library
find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(tqstests PUBLIC ${MATH_LIBRARY})
endif()

# Interprocedural optimizations
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if (result)
        set_target_properties(tqstests PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

# GNU Scientific Library
find_package(GSL REQUIRED)
if (GSL_FOUND)
    if (GSL_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "tqs needs GSL >= 2.0")
    endif()

    target_include_directories(tqstests PUBLIC "${GSL_INCLUDE_DIRS}")
    target_link_libraries(tqstests PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(tqstests PUBLIC "${GSL_CBLAS_LIBRARY}")
endif()

# SOFTLib
target_link_libraries(tqstests PUBLIC softlib)
#message(INFO "Looking for softlib...")
#find_package(SOFTLIB REQUIRED)
#if (SOFTLIB_FOUND)
#    target_include_directories(tqs PUBLIC "${SOFTLIB_INCLUDE_DIRS}")
#    target_link_libraries(tqs PUBLIC "${SOFTLIB_LIBRARIES}")
#endif()

# PETSc
find_package(PETSc COMPONENTS CXX REQUIRED)
if (PETSC_FOUND)
    target_include_directories(tqstests PUBLIC "${PETSC_INCLUDES}")
    target_link_libraries(tqstests PUBLIC "${PETSC_LIBRARIES}")
    add_definitions(${PETSC_DEFINITIONS})
endif()

