
set(dreami_files
    "${PROJECT_SOURCE_DIR}/iface/Main.cpp"
)

add_executable(dreami ${dreami_files})
target_link_libraries(dreami PUBLIC dream)
set_target_properties(dreami PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

