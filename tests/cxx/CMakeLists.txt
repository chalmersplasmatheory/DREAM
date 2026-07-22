
set(dreamtests_core
    "${PROJECT_SOURCE_DIR}/tests/cxx/runtest.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/UnitTest.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/UnitTest.grids.cpp"
)

set(dreamtests_dream
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/DREAM/AvalancheSourceRP.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/DREAM/BoundaryFlux.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/DREAM/IonRateEquation.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/DREAM/MeanExcitationEnergy.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/DREAM/RunawayFluid.cpp"
)

set(dreamtests_fvm
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/AdvectionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/AdvectionDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/AnalyticBRadialGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/DiffusionTerm.cpp"    
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/EquationTerm.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/GeneralAdvectionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/GeneralAdvectionDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/GeneralDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/Grid.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/Interpolator1D.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/Interpolator3D.cpp"
    "${PROJECT_SOURCE_DIR}/tests/cxx/tests/FVM/PXiExternalKineticKinetic.cpp"
)

add_executable(dreamtests ${dreamtests_core} ${dreamtests_dream} ${dreamtests_fvm})
target_include_directories(dreamtests PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/tests/cxx/include")
target_link_libraries(dreamtests PUBLIC fvm)
target_link_libraries(dreamtests PUBLIC dream)

# Require C++17
set_target_properties(dreamtests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

######################
# MODULES
######################
# Math library
find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(dreamtests PUBLIC ${MATH_LIBRARY})
endif()

# Interprocedural optimizations
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if (result)
        set_target_properties(dreamtests PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

# GNU Scientific Library
find_package(GSL REQUIRED)
if (GSL_FOUND)
    if (GSL_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "dream needs GSL >= 2.0")
    endif()

    target_include_directories(dreamtests PUBLIC "${GSL_INCLUDE_DIRS}")
    target_link_libraries(dreamtests PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(dreamtests PUBLIC "${GSL_CBLAS_LIBRARY}")
endif()

# SOFTLib
target_link_libraries(dreamtests PUBLIC softlib)
#message(INFO "Looking for softlib...")
#find_package(SOFTLIB REQUIRED)
#if (SOFTLIB_FOUND)
#    target_include_directories(dream PUBLIC "${SOFTLIB_INCLUDE_DIRS}")
#    target_link_libraries(dream PUBLIC "${SOFTLIB_LIBRARIES}")
#endif()

# PETSc
find_package(PETSc COMPONENTS CXX REQUIRED)
if (PETSC_FOUND)
    target_include_directories(dreamtests PUBLIC "${PETSC_INCLUDES}")
    target_link_libraries(dreamtests PUBLIC "${PETSC_LIBRARIES}")
    add_definitions(${PETSC_DEFINITIONS})
endif()

