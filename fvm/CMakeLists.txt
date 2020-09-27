
set(fvm_core
    "${PROJECT_SOURCE_DIR}/fvm/BlockMatrix.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/DurationTimer.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Init.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Interpolator1D.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Interpolator3D.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Matrix.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/NormEvaluator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/UnknownQuantity.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/UnknownQuantityHandler.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/QuantityData.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Solvers/MILU.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Solvers/MIKSP.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/TimeKeeper.cpp"
)
set(fvm_core_headers
    "${PROJECT_SOURCE_DIR}/include/FVM/BlockMatrix.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/FVMException.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Matrix.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/MatrixInverter.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Solvers/MILU.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Solvers/MIKSP.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/fluxGridType.enum.hpp"
)
set(fvm_equation
    "${PROJECT_SOURCE_DIR}/fvm/Equation/AdvectionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/AdvectionDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/AdvectionInterpolationCoefficient.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/BoundaryConditions/PXiExternalKineticKinetic.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/BoundaryConditions/PXiExternalLoss.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/BoundaryConditions/PInternalBoundaryCondition.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/BoundaryConditions/XiInternalBoundaryCondition.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/ConstantParameter.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/DiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/Operator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/EquationTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/DiagonalTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/DiagonalLinearTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/DiagonalQuadraticTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/DiagonalComplexTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/ScalarLinearTerm.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/MomentQuantity.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/PredeterminedParameter.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/PrescribedParameter.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Equation/LinearTransientTerm.cpp"
)
set(fvm_equation_headers
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/Operator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/EquationTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/AdvectionTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/AdvectionDiffusionTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/DiffusionTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/TransientTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/IdentityTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Equation/TransientTerm.hpp"
)
set(fvm_grid
    "${PROJECT_SOURCE_DIR}/fvm/Grid/CylindricalRadialGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/AnalyticBRadialGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/EmptyMomentumGrid.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/EmptyRadialGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/Grid.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/MomentumGrid.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/PUniformGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/PBiUniformGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/PXiMomentumGrid.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/PXiMomentumGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/PXiGrid/XiUniformGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/RadialGrid.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/RadialGridGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/FluxSurfaceAverager.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/FluxSurfaceQuantity.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/BounceAverager.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/BounceSurfaceQuantity.cpp"
    "${PROJECT_SOURCE_DIR}/fvm/Grid/BounceSurfaceMetric.cpp"
)
set(fvm_grid_headers
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/CylindricalRadialGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/Grid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/EmptyMomentumGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/EmptyRadialGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/MomentumGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/MomentumGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PUniformGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PBiUniformGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PXiMomentumGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/PXiMomentumGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/PXiGrid/XiUniformGridGenerator.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/RadialGrid.hpp"
    "${PROJECT_SOURCE_DIR}/include/FVM/Grid/RadialGridGenerator.hpp"
)

add_library(fvm STATIC ${fvm_core} ${fvm_equation} ${fvm_grid} ${fvm_core_headers} ${fvm_equation_headers} ${fvm_grid_headers})
target_include_directories(fvm PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/extern/softlib/include")

set_target_properties(fvm PROPERTIES
    # Require C++17
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if (POS_INDP_CODE)
    # Make code position-independent (required for linking to Python interface)
    set_property(TARGET fvm PROPERTY POSITION_INDEPENDENT_CODE YES)
endif ()

# Link with softlib
target_link_libraries(fvm PUBLIC softlib)

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
        message(FATAL_ERROR "DREAM needs GSL >= 2.0")
    endif()

    target_include_directories(fvm PUBLIC "${GSL_INCLUDE_DIRS}")
    target_link_libraries(fvm PUBLIC "${GSL_LIBRARIES}")
    target_link_libraries(fvm PUBLIC "${GSL_CBLAS_LIBRARY}")
endif()

# MPI (just as a dependency for PETSc)
find_package(MPI COMPONENTS CXX)
if (MPI_CXX_FOUND)
    target_include_directories(fvm PUBLIC ${MPI_CXX_INCLUDE_DIRS})
	set_target_properties(fvm PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_OPTIONS} ${MPI_CXX_COMPILE_DEFINITIONS}")
	set_target_properties(fvm PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
	target_link_libraries(fvm PUBLIC "${MPI_CXX_LIBRARIES}")
endif (MPI_CXX_FOUND)

# PETSc
find_package(PETSc COMPONENTS CXX REQUIRED)
if (PETSC_FOUND)
    target_include_directories(fvm PUBLIC "${PETSC_INCLUDES}")
    target_link_libraries(fvm PUBLIC "${PETSC_LIBRARIES}")
    add_definitions(${PETSC_DEFINITIONS})
endif()

