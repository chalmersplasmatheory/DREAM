
set(tqs_core
    "${PROJECT_SOURCE_DIR}/src/Init.cpp"
)
set(tqs_headers
    "${PROJECT_BINARY_DIR}/include/config.h"
    "${PROJECT_SOURCE_DIR}/include/TQS/Init.h"
)

add_library(tqs STATIC ${tqs_core} ${tqs_headers})
target_include_directories(tqs PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include")

# Require C++17
set_target_properties(tqs PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

######################
# MODULES
######################
# Math library
find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(tqs PUBLIC ${MATH_LIBRARY})
endif()

# Interprocedural optimizations
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if (result)
        set_target_properties(tqs PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

# GNU Scientific Library
find_package(GSL REQUIRED)
if (GSL_FOUND)
    if (GSL_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "tqs needs GSL >= 2.0")
    endif()

    target_include_directories(tqs PUBLIC "${GSL_INCLUDE_DIRS}")
    target_link_libraries(tqs PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(tqs PUBLIC "${GSL_CBLAS_LIBRARY}")
endif()

find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
    target_include_directories(tqs PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(tqs PUBLIC ${HDF5_LIBRARIES})
else (HDF5_FOUND)
    message(FATAL_ERROR "No HDF5 C++ library was found on this system")
endif (HDF5_FOUND)

# SOFTLib
target_link_libraries(tqs PUBLIC softlib)
#message(INFO "Looking for softlib...")
#find_package(SOFTLIB REQUIRED)
#if (SOFTLIB_FOUND)
#    target_include_directories(tqs PUBLIC "${SOFTLIB_INCLUDE_DIRS}")
#    target_link_libraries(tqs PUBLIC "${SOFTLIB_LIBRARIES}")
#endif()

# PETSc
find_package(PETSc COMPONENTS CXX REQUIRED)
if (PETSC_FOUND)
    target_include_directories(tqs PUBLIC "${PETSC_INCLUDES}")
    target_link_libraries(tqs PUBLIC "${PETSC_LIBRARIES}")
    add_definitions(${PETSC_DEFINITIONS})
endif()

