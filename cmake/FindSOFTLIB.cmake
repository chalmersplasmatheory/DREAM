#.rst:
# FindSOFTLib
# -----------
#
# Finds SOFTLib, both binaries and headers.

set(SOFTLIB_INCLUDE_RELPATH "extern/softlib/include/")
set(SOFTLIB_LIBRARY_RELPATH "extern/softlib/build/src/")

find_path(SOFTLIB_INCLUDE_DIR softlib/config.h.in
    HINTS ${SOFTLIB_INCLUDE_RELPATH})
find_library(SOFTLIB_LIBRARY softlib libsoftlib
    HINTS ${SOFTLIB_LIBRARY_RELPATH})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SOFTLIB DEFAULT_MSG SOFTLIB_LIBRARY SOFTLIB_INCLUDE_DIR)
mark_as_advanced(SOFTLIB_INCLUDE_DIR SOFTLIB_LIBRARY)

set(SOFTLIB_LIBRARIES ${SOFTLIB_LIBRARY})
set(SOFTLIB_INCLUDE_DIRS ${SOFTLIB_INCLUDE_DIR})
