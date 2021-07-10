
set(dreampy_files
    "${PROJECT_SOURCE_DIR}/pyface/dreampy.cpp"
    "${PROJECT_SOURCE_DIR}/pyface/numpy.cpp"
    "${PROJECT_SOURCE_DIR}/pyface/settings.cpp"
    "${PROJECT_SOURCE_DIR}/pyface/SFile_Python.cpp"
)

add_library(dreampy SHARED ${dreampy_files})
target_link_libraries(dreampy PUBLIC dream)
set_target_properties(dreampy PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Find the Python API and link to it
find_package(Python3 REQUIRED COMPONENTS Development NumPy)
if (Python3_FOUND)
    message(STATUS "${Python3_NumPy_INCLUDE_DIRS}")
    target_include_directories(dreampy PUBLIC "${Python3_INCLUDE_DIRS}" "${Python3_NumPy_INCLUDE_DIRS}")
    target_link_libraries(dreampy PUBLIC "${Python3_LIBRARIES}")
endif ()

