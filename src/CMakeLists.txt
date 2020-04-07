
set(dream_core
    "${PROJECT_SOURCE_DIR}/src/Constants.cpp"
    "${PROJECT_SOURCE_DIR}/src/EquationSystem.cpp"
    "${PROJECT_SOURCE_DIR}/src/IO.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init.cpp"
    "${PROJECT_SOURCE_DIR}/src/Simulation.cpp"
    "${PROJECT_SOURCE_DIR}/src/UnknownQuantityEquation.cpp"
)
set(dream_equations
    "${PROJECT_SOURCE_DIR}/src/Equations/CollisionQuantityHandler.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/SlowingDownTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/ElectricFieldTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/ElectricFieldDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/PitchScatterTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/EnergyDiffusionTerm.cpp"
)
set(dream_settings
    "${PROJECT_SOURCE_DIR}/src/Settings/Constants.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/DefineOptions.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/n_cold.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/EquationSystem.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/KineticGrid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Process.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/RadialGrid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Settings.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/SFile.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Solver.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/TimeStepper.cpp"
)
set(dream_headers
    "${PROJECT_SOURCE_DIR}/include/DREAM/IO.hpp"
    "${PROJECT_SOURCE_DIR}/include/DREAM/Init.h"
    "${PROJECT_SOURCE_DIR}/include/DREAM/Settings/Settings.hpp"
    "${PROJECT_SOURCE_DIR}/include/DREAM/Settings/SFile.hpp"
)

set(dream_solvers
    "${PROJECT_SOURCE_DIR}/src/Solver/Solver.cpp"
    "${PROJECT_SOURCE_DIR}/src/Solver/SolverLinearlyImplicit.cpp"
    "${PROJECT_SOURCE_DIR}/src/Solver/SolverSNES.cpp"
    "${PROJECT_SOURCE_DIR}/src/Solver/SNES_helpers.cpp"
)

add_library(dream STATIC ${dream_core} ${dream_settings} ${dream_solvers} ${dream_equations} ${dream_headers})
target_include_directories(dream PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(dream PUBLIC fvm)

# Require C++17
set_target_properties(dream PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

######################
# MODULES
######################
# Math library
find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(dream PUBLIC ${MATH_LIBRARY})
endif()

# Interprocedural optimizations
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if (result)
        set_target_properties(dream PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

# GNU Scientific Library
find_package(GSL REQUIRED)
if (GSL_FOUND)
    if (GSL_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "dream needs GSL >= 2.0")
    endif()

    target_include_directories(dream PUBLIC "${GSL_INCLUDE_DIRS}")
    target_link_libraries(dream PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(dream PUBLIC "${GSL_CBLAS_LIBRARY}")
endif()

find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
    target_include_directories(dream PUBLIC ${HDF5_INCLUDE_DIRS})
    target_link_libraries(dream PUBLIC ${HDF5_LIBRARIES})
else (HDF5_FOUND)
    message(FATAL_ERROR "No HDF5 C++ library was found on this system")
endif (HDF5_FOUND)

# SOFTLib
target_link_libraries(dream PUBLIC softlib)
#message(INFO "Looking for softlib...")
#find_package(SOFTLIB REQUIRED)
#if (SOFTLIB_FOUND)
#    target_include_directories(dream PUBLIC "${SOFTLIB_INCLUDE_DIRS}")
#    target_link_libraries(dream PUBLIC "${SOFTLIB_LIBRARIES}")
#endif()

# PETSc
find_package(PETSc COMPONENTS CXX REQUIRED)
if (PETSC_FOUND)
    target_include_directories(dream PUBLIC "${PETSC_INCLUDES}")
    target_link_libraries(dream PUBLIC "${PETSC_LIBRARIES}")
    add_definitions(${PETSC_DEFINITIONS})
endif()

