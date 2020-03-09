
set(tqsi_files
    "${PROJECT_SOURCE_DIR}/iface/Main.cpp"
)

add_executable(tqsi ${tqsi_files})
target_link_libraries(tqsi PRIVATE tqs)
set_target_properties(tqsi PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

