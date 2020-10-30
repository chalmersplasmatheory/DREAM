
set(dream_core
    "${PROJECT_SOURCE_DIR}/src/ADAS.cpp"
    "${PROJECT_SOURCE_DIR}/src/ADASRateInterpolator.cpp"
    "${PROJECT_SOURCE_DIR}/src/Atomics/adasdata.cpp"
    "${PROJECT_SOURCE_DIR}/src/Atomics/nistdata_binding.cpp"
    "${PROJECT_SOURCE_DIR}/src/Atomics/nistdata_ionization.cpp"
    "${PROJECT_SOURCE_DIR}/src/Constants.cpp"
    "${PROJECT_SOURCE_DIR}/src/ConvergenceChecker.cpp"
    "${PROJECT_SOURCE_DIR}/src/EquationSystem/EquationSystem.cpp"
    "${PROJECT_SOURCE_DIR}/src/EquationSystem/Info.cpp"
    "${PROJECT_SOURCE_DIR}/src/EquationSystem/Save.cpp"
    "${PROJECT_SOURCE_DIR}/src/EqsysInitializer.cpp"
    "${PROJECT_SOURCE_DIR}/src/IO.cpp"
    "${PROJECT_SOURCE_DIR}/src/IonHandler.cpp"
    "${PROJECT_SOURCE_DIR}/src/IonInterpolator1D.cpp"
    "${PROJECT_SOURCE_DIR}/src/Init.cpp"
    "${PROJECT_SOURCE_DIR}/src/NIST.cpp"
    "${PROJECT_SOURCE_DIR}/src/OutputGenerator.cpp"
    "${PROJECT_SOURCE_DIR}/src/OutputGeneratorSFile.cpp"
    "${PROJECT_SOURCE_DIR}/src/OtherQuantityHandler.cpp"
    "${PROJECT_SOURCE_DIR}/src/PostProcessor.cpp"
    "${PROJECT_SOURCE_DIR}/src/Simulation.cpp"
    "${PROJECT_SOURCE_DIR}/src/TimeStepper/TimeStepperAdaptive.cpp"
    "${PROJECT_SOURCE_DIR}/src/TimeStepper/TimeStepperConstant.cpp"
    "${PROJECT_SOURCE_DIR}/src/UnknownQuantityEquation.cpp"
)
set(dream_equations
    "${PROJECT_SOURCE_DIR}/src/Equations/AnalyticDistributionRE.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/CollisionQuantityHandler.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/CollisionQuantity.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/CollisionFrequency.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/ConnorHastie.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/CoulombLogarithm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/DreicerNeuralNetwork.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/EffectiveCriticalField.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/FluidSourceTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/BindingEnergyTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/CollisionalEnergyTransferKineticTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/ComptonRateTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/CurrentDensityFromDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/DensityFromBoundaryFluxPXI.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/DensityFromDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/DreicerRateTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/FreeElectronDensityTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/HeatTransportRechesterRosenbluth.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/HotTailCurrentDensityFromDistributionFunction.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/IonKineticIonizationTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/IonPrescribedParameter.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/IonRateEquation.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/IonTransientTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/FreeElectronDensityTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/TotalElectronDensityTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/RadiatedPowerTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Fluid/TotalElectronDensityTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/AvalancheSourceRP.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/BCIsotropicSourcePXi.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/ElectricFieldTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/ElectricFieldDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/EnergyDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/PitchScatterTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/ParticleSourceTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/RechesterRosenbluthTransport.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/SlowingDownTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/Kinetic/SynchrotronTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/ParallelDiffusionFrequency.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/PitchScatterFrequency.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/PoloidalFlux/AmperesLawBoundaryAtRMax.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/PoloidalFlux/HyperresistiveDiffusionTerm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/SlowingDownFrequency.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/RunawayFluid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/TransportPrescribed.cpp"
    "${PROJECT_SOURCE_DIR}/src/Equations/TransportBC.cpp"
)
set(dream_equation_headers
    "${PROJECT_SOURCE_DIR}/include/DREAM/Equations/Fluid/CurrentFromConductivityTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/DREAM/Equations/Fluid/AvalancheGrowthTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/DREAM/Equations/Fluid/PredictedOhmicCurrentFromDistributionTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/DREAM/Equations/Fluid/OhmicHeatingTerm.hpp"
    "${PROJECT_SOURCE_DIR}/include/DREAM/Equations/Scalar/WallCurrentTerms.hpp"
)
set(dream_settings
    "${PROJECT_SOURCE_DIR}/src/Settings/CollisionQuantityHandler.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/RunawayFluid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Constants.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/DefineOptions.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/distribution.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/E_field.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/f_hot.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/f_re.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/ions.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/n_cold.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/n_hot.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/j_hot.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/j_re.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/j_ohm.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/j_tot.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/psi_p.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/n_re.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/n_tot.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Equations/T_cold.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/EquationSystem.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/KineticGrid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/LoadPrescribedData.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/OtherQuantityHandler.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Process.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/RadialGrid.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Settings.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/SFile.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Solver.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/TimeStepper.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Transport.cpp"
    "${PROJECT_SOURCE_DIR}/src/Settings/Tolerances.cpp"
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
    "${PROJECT_SOURCE_DIR}/src/Solver/SolverNonLinear.cpp"
    "${PROJECT_SOURCE_DIR}/src/Solver/NumericalJacobian.cpp"
)

add_library(dream STATIC ${dream_core} ${dream_settings} ${dream_solvers} ${dream_equations} ${dream_equation_headers} ${dream_headers})
target_include_directories(dream PUBLIC "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(dream PUBLIC fvm)

set_target_properties(dream PROPERTIES
    # Require C++17
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)

if (POS_INDP_CODE)
    # Make code position-independent (required for linking to Python interface)
    set_property(TARGET dream PROPERTY POSITION_INDEPENDENT_CODE YES)
endif ()

# Command for generating ADAS data (requires Python 3)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
if (Python3_FOUND)
    add_custom_command(OUTPUT "${PROJECT_SOURCE_DIR}/src/Atomics/adasdata.cpp" COMMAND "${Python3_EXECUTABLE}" "${PROJECT_SOURCE_DIR}/tools/get_adas.py" WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/tools" DEPENDS "${PROJECT_SOURCE_DIR}/tools/elements.json")
endif ()

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

# MPI (just as a dependency for PETSc)
find_package(MPI COMPONENTS CXX)
if (MPI_CXX_FOUND)
    target_include_directories(dream PUBLIC ${MPI_CXX_INCLUDE_DIRS})
	set_target_properties(dream PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_OPTIONS} ${MPI_CXX_COMPILE_DEFINITIONS}")
	set_target_properties(dream PROPERTIES LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
	target_link_libraries(dream PUBLIC "${MPI_CXX_LIBRARIES}")
endif (MPI_CXX_FOUND)

# PETSc
find_package(PETSc COMPONENTS CXX REQUIRED)
if (PETSC_FOUND)
    target_include_directories(dream PUBLIC "${PETSC_INCLUDES}")
    target_link_libraries(dream PUBLIC "${PETSC_LIBRARIES}")
    add_definitions(${PETSC_DEFINITIONS})
endif()

