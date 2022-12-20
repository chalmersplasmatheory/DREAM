# - Determines if PETSc was compiled with MPI dependencies.
#
# Written by: Mathias Hoppe, December 2022
#

if (__check_petsc_requires_mpi)
	return ()
endif()
set(__check_petsc_requires_mpi YES)

function(check_petsc_requires_mpi _petsc_mpi_var)
	set(_incdir "${PETSC_DIR}/${PETSC_ARCH}/include")
	set(_incfile "${_incdir}/petscconfiginfo.h")
	if (EXISTS "${_incfile}")
		file(READ "${_incfile}" PETSC_CONFIG_INFO)

		if (PETSC_CONFIG_INFO MATCHES "--with-mpi=0")
			set(${_petsc_mpi_var} OFF PARENT_SCOPE)
			return()
		else ()
			set(${_petsc_mpi_var} ON PARENT_SCOPE)
			return()
		endif ()
	endif()
endfunction()

