
set(fvm_core
    "${PROJECT_SOURCE_DIR}/fvm/EquationSystem.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Matrix.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Solvers/MILU.cpp"
)
set(fvm_core_headers
    "${PROJECT_BINARY_DIR}/include/config.h"
    "${PROJECT_SOURCE_DIR}/include/FVM/EquationSystem.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/FVMException.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Matrix.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/MatrixInverter.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Solvers/MILU.hpp"
)
set(fvm_equation
    "${PROJECT_SOURCE_DIR}/fvm/Equation/Equation.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/EquationTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/AdvectionTerm.cpp"
)
set(fvm_equation_headers
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/Equation.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/EquationTerm.hpp"
)
set(fvm_grid
    "${PROJECT_SOURCE_DIR}/fvm/Grid/CylindricalRadialGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/MomentumGrid.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/PUniformGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/PXiMomentumGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/XiUniformGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/RadialGrid.cpp"
)
set(fvm_grid_headers
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/CylindricalRadialGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/MomentumGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/MomentumGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/RadialGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/RadialGridGenerator.hpp"
)

add_library(fvm STATIC ${fvm_core} ${fvm_equation} ${fvm_grid} ${fvm_core_headers} ${fvm_equation_headers} ${fvm_grid_headers})
target_include_directories(fvm PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/extern/softlib/include")

# Require C++17
set_target_properties(fvm PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

######################
# MODULES
######################
# Math library
find_library(MATH_LIBRARY m)
if (MATH_LIBRARY)
    target_link_libraries(fvm PUBLIC ${MATH_LIBRARY})
endif()

# Interprocedural optimizations
if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if (result)
        set_target_properties(fvm PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

# GNU Scientific Library
find_package(GSL REQUIRED)
if (GSL_FOUND)
    if (GSL_VERSION VERSION_LESS 2.0)
        message(FATAL_ERROR "TQS needs GSL >= 2.0")
    endif()

    target_include_directories(fvm PUBLIC "${GSL_INCLUDE_DIRS}")
    target_link_libraries(fvm PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(fvm PUBLIC "${GSL_CBLAS_LIBRARY}")
endif()

# PETSc
find_package(PETSc COMPONENTS CXX REQUIRED)
if (PETSC_FOUND)
    target_include_directories(fvm PUBLIC "${PETSC_INCLUDES}")
    target_link_libraries(fvm PUBLIC "${PETSC_LIBRARIES}")
    add_definitions(${PETSC_DEFINITIONS})
endif()

