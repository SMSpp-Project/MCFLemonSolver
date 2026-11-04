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
#    and, after it, from the environment variable LEMON_ROOT and from the     #
#    directories where INSTALL.sh and Homebrew put abLEMON, the maintained    #
#    fork of LEMON, before the directories of the system.                     #
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

# Drop the cached paths if one of them does not exist any more, so that
# they are all searched again from the same place
set(_LEMON_cached LEMON_INCLUDE_DIR LEMON_LIBRARY LEMON_LIBRARY_DEBUG
        LEMON_VERSION)
foreach (_var ${_LEMON_cached})
    foreach (_path ${${_var}})
        if (IS_ABSOLUTE "${_path}" AND NOT EXISTS "${_path}")
            set(_LEMON_stale TRUE)
        endif ()
    endforeach ()
endforeach ()
# the hints: LEMON_ROOT, the environment, and where abLEMON is installed by
# INSTALL.sh (/opt/abLEMON) and by Homebrew
set(_LEMON_hints ${LEMON_ROOT} $ENV{LEMON_ROOT} /opt/abLEMON
        /opt/homebrew/opt/ablemon /usr/local/opt/ablemon)

# A LEMON found before is searched again if it is not under the LEMON_ROOT
# given now, or if it is not abLEMON (lemon/bits/mcf_certificate.h) while
# one of the hints holds abLEMON, as when the LEMON of the system was found
# before abLEMON was installed
if (LEMON_INCLUDE_DIR)
    if (LEMON_ROOT)
        file(REAL_PATH "${LEMON_ROOT}" _LEMON_root)
        file(REAL_PATH "${LEMON_INCLUDE_DIR}" _LEMON_inc)
        string(FIND "${_LEMON_inc}" "${_LEMON_root}" _LEMON_at)
        if (NOT _LEMON_at EQUAL 0)
            set(_LEMON_stale TRUE)
        endif ()
    endif ()
    if (NOT EXISTS "${LEMON_INCLUDE_DIR}/lemon/bits/mcf_certificate.h")
        foreach (_LEMON_hint ${_LEMON_hints})
            if (EXISTS "${_LEMON_hint}/include/lemon/bits/mcf_certificate.h")
                set(_LEMON_stale TRUE)
            endif ()
        endforeach ()
    endif ()
endif ()

if (_LEMON_stale)
    message(STATUS "Cached LEMON paths not found, searching again")
    foreach (_var ${_LEMON_cached})
        unset(${_var} CACHE)
    endforeach ()
endif ()
unset(_LEMON_cached)
unset(_LEMON_stale)
unset(_var)
unset(_path)

# The cached version is the one of the include directory it was read from
if (NOT "${LEMON_INCLUDE_DIR}" STREQUAL
        "${_LEMON_VERSION_INCLUDE_DIR}")
    unset(LEMON_VERSION CACHE)
endif ()

# Check if already in cache
if (LEMON_INCLUDE_DIR AND LEMON_LIBRARY AND LEMON_LIBRARY_DEBUG
        AND LEMON_VERSION)
    set(LEMON_FOUND TRUE)
endif ()

if (NOT LEMON_FOUND)

    # ----- Find the LEMON include directory -------------------------------- #
    # the hints are searched before the system directories, so that the
    # LEMON they point to wins over another one installed in the system
    find_path(LEMON_INCLUDE_DIR
            NAMES lemon/config.h lemon/network_simplex.h
            HINTS ${_LEMON_hints}
            PATH_SUFFIXES include
            DOC "LEMON include directory.")

    # ----- Find the LEMON library ------------------------------------------ #
    # the library is named lemon by the distributions and emon by the CMake
    # build of LEMON: NAMES_PER_DIR looks for both names in a directory
    # before the next one, so that the emon of LEMON_ROOT is not passed over
    # for the lemon of the system
    if (UNIX)
        # the library of the directory where the headers were found first
        get_filename_component(_LEMON_lib_hints "${LEMON_INCLUDE_DIR}"
                DIRECTORY)
        set(_LEMON_lib_hints ${_LEMON_lib_hints}/lib)
        foreach (_LEMON_hint ${_LEMON_hints})
            list(APPEND _LEMON_lib_hints ${_LEMON_hint}/lib)
        endforeach ()
        find_library(LEMON_LIBRARY
                NAMES lemon emon
                NAMES_PER_DIR
                HINTS ${_LEMON_lib_hints}
                DOC "LEMON library.")

        set(LEMON_LIBRARY_DEBUG ${LEMON_LIBRARY}
                CACHE FILEPATH "LEMON debug library." FORCE)
    elseif (WIN32)
        find_library(LEMON_LIBRARY
                NAMES lemon emon
                NAMES_PER_DIR
                PATHS
                ${LEMON_ROOT}/lib
                ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/lib
                $ENV{LIBRARY_LIB}
                NO_DEFAULT_PATH
                DOC "LEMON library.")

        find_library(LEMON_LIBRARY_DEBUG
                NAMES lemon emon
                NAMES_PER_DIR
                PATHS
                ${LEMON_ROOT}/debug/lib
                ${VCPKG_INSTALLED_DIR}/${VCPKG_TARGET_TRIPLET}/debug/lib
                NO_DEFAULT_PATH
                DOC "LEMON debug library.")

        # Release-only distributions (e.g. conda-forge) ship no debug build:
        # fall back to the release library so a Release configure succeeds.
        if (NOT LEMON_LIBRARY_DEBUG)
            set(LEMON_LIBRARY_DEBUG ${LEMON_LIBRARY}
                    CACHE FILEPATH "LEMON debug library." FORCE)
        endif ()
    endif ()

    # ----- Parse the version ----------------------------------------------- #
    # lemon/config.h carries a single string macro, e.g.,
    # #define LEMON_VERSION "1.3.1".
    if (LEMON_INCLUDE_DIR AND EXISTS "${LEMON_INCLUDE_DIR}/lemon/config.h")
        file(STRINGS
                "${LEMON_INCLUDE_DIR}/lemon/config.h"
                _LEMON_version_line
                REGEX "#define[ \t]+LEMON_VERSION[ \t]+\"[^\"]")

        if (_LEMON_version_line)
            string(REGEX REPLACE ".*\"([^\"]+)\".*" "\\1"
                    LEMON_VERSION "${_LEMON_version_line}")
        endif ()
        unset(_LEMON_version_line)
    endif ()

    # ----- Handle the standard arguments ----------------------------------- #
    # The following macro manages the QUIET, REQUIRED and version-related
    # options passed to find_package(). It also sets <PackageName>_FOUND if
    # REQUIRED_VARS are set.
    # REQUIRED_VARS should be cache entries and not output variables. See:
    # https://cmake.org/cmake/help/latest/module/
    #   FindPackageHandleStandardArgs.html
    # The debug library is optional (it falls back to the release one above),
    # so it is deliberately kept out of REQUIRED_VARS.
    find_package_handle_standard_args(
            LEMON
            REQUIRED_VARS LEMON_LIBRARY LEMON_INCLUDE_DIR
            VERSION_VAR LEMON_VERSION)
endif ()

# A found package keeps its version in the cache, with the include
# directory it was read from
if (LEMON_FOUND)
    set(LEMON_VERSION "${LEMON_VERSION}" CACHE INTERNAL
            "LEMON version.")
    set(_LEMON_VERSION_INCLUDE_DIR "${LEMON_INCLUDE_DIR}"
            CACHE INTERNAL "Include directory of the LEMON version.")
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
                IMPORTED_LOCATION_DEBUG "${LEMON_LIBRARY_DEBUG}"
                INTERFACE_INCLUDE_DIRECTORIES "${LEMON_INCLUDE_DIRS}"
                INTERFACE_LINK_LIBRARIES "${CMAKE_THREAD_LIBS_INIT}")
    endif ()
endif ()

# Variables marked as advanced are not displayed in CMake GUIs, see:
# https://cmake.org/cmake/help/latest/command/mark_as_advanced.html
mark_as_advanced(LEMON_INCLUDE_DIR
        LEMON_LIBRARY
        LEMON_LIBRARY_DEBUG
        LEMON_VERSION)

# --------------------------------------------------------------------------- #
