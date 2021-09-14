
set(dreami_files
    "${PROJECT_SOURCE_DIR}/iface/Main.cpp"
)

add_executable(dreami ${dreami_files})
target_link_libraries(dreami PUBLIC dream)
set_target_properties(dreami PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    include(CheckIPOSupported)
    # Clang does not support IPO (although Clang++ does), so we need
    # to explicitly check both C and C++ (since Clang is used for linking)
    check_ipo_supported(RESULT result LANGUAGES C CXX)
    if (result)
        set_target_properties(dreami PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif ()

