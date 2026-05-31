# --------------------------------------------------------------------------- #
#    CMake find module for LEMON                                              #
#                                                                             #
#    This module finds LEMON include directories and libraries.               #
#    Use it by invoking find_package() with the form:                         #
#                                                                             #
#        find_package(LEMON [version] [EXACT] [REQUIRED])                     #
#                                                                             #
#    The results are stored in the following variables:                       #
#                                                                             #
#        LEMON_FOUND         - True if headers are found                      #
#        LEMON_INCLUDE_DIRS  - Include directories                            #
#        LEMON_LIBRARIES     - Libraries to be linked                         #
#        LEMON_VERSION       - Version number                                 #
#                                                                             #
#    This module reads hints about search locations from variables:           #
#                                                                             #
#        LEMON_ROOT          - Custom path to LEMON                           #
#                                                                             #
#    The following IMPORTED target is also defined:                           #
#                                                                             #
#        LEMON::LEMON                                                         #
#                                                                             #
#    This find module is provided because LEMON ships its CMake config under  #
#    an inconsistent, distributor-dependent name and exposes the library      #
#    under more than one output name, so find_package(... CONFIG) is not      #
#    portable.                                                                #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_THREAD_LIBS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
# LEMON's lock primitives (lemon/bits/lock.h) may pull in pthreads.
find_package(Threads QUIET)

# Check if already in cache
if (LEMON_INCLUDE_DIR AND LEMON_LIBRARY AND LEMON_VERSION)
    set(LEMON_FOUND TRUE)
endif ()

if (NOT LEMON_FOUND)

    # ----- Find the LEMON include directory -------------------------------- #
    find_path(LEMON_INCLUDE_DIR
            NAMES lemon/config.h lemon/network_simplex.h
            PATHS ${LEMON_ROOT}
            PATH_SUFFIXES include
            DOC "LEMON include directory.")

    # ----- Find the LEMON library ------------------------------------------ #
    # The library is named "lemon" on Debian (liblemon) but "emon" upstream
    # (libemon, the historical COIN-OR output name): look for both.
    find_library(LEMON_LIBRARY
            NAMES lemon emon
            PATHS ${LEMON_ROOT}/lib
            DOC "LEMON library.")

    # ----- Parse the version ----------------------------------------------- #
    # lemon/config.h carries a single string macro: #define LEMON_VERSION "1.3.1"
    if (LEMON_INCLUDE_DIR AND EXISTS "${LEMON_INCLUDE_DIR}/lemon/config.h")
        file(STRINGS
                "${LEMON_INCLUDE_DIR}/lemon/config.h"
                _LEMON_version_line REGEX "#define[ \t]+LEMON_VERSION[ \t]+\"")

        string(REGEX REPLACE ".*#define[ \t]+LEMON_VERSION[ \t]+\"([^\"]+)\".*" "\\1"
                LEMON_VERSION "${_LEMON_version_line}")
        unset(_LEMON_version_line)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/FindPackageHandleStandardArgs.html
    find_package_handle_standard_args(
            LEMON
            REQUIRED_VARS LEMON_LIBRARY LEMON_INCLUDE_DIR
            VERSION_VAR LEMON_VERSION)
endif ()

# ----- Export the target --------------------------------------------------- #
if (LEMON_FOUND)
    set(LEMON_INCLUDE_DIRS ${LEMON_INCLUDE_DIR})
    set(LEMON_LIBRARIES ${LEMON_LIBRARY})

    if (CMAKE_THREAD_LIBS_INIT)
        set(LEMON_LIBRARIES ${LEMON_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
    endif ()

    if (NOT TARGET LEMON::LEMON)
        add_library(LEMON::LEMON UNKNOWN IMPORTED)
        set_target_properties(
                LEMON::LEMON PROPERTIES
                IMPORTED_LOCATION "${LEMON_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${LEMON_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(LEMON_INCLUDE_DIR
        LEMON_LIBRARY
        LEMON_VERSION)

# --------------------------------------------------------------------------- #
