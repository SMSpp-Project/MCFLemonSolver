# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ MCFLemonSolver                                 #
#                                                                             #
#    This file allows one to build the library using CMake.                   #
#    To do so, you can use the following commands:                            #
#                                                                             #
#        $ cmake -S <source-path> -B <build-path>                             #
#        $ cmake --build <build-path>                                         #
#                                                                             #
#    The following command also installs the library in the system:           #
#                                                                             #
#        $ cmake --build <build-path> --target install                        #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #
cmake_minimum_required(VERSION 3.21)
cmake_policy(VERSION 3.15)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DeriveVersion.cmake)
smspp_derive_version(SMSPP_MODULE_VERSION)

project(MCFLemonSolver
        VERSION ${SMSPP_MODULE_VERSION}
        DESCRIPTION "SMS++ MCFLemonSolver"
        HOMEPAGE_URL https://gitlab.com/smspp/mcflemonsolver
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})

# C-safe form of modName for the preprocessor macros in <modName>Config.h
string(MAKE_C_IDENTIFIER "${modName}" modMacro)
set(modNamespace "SMS++")

# Find out if it's being called by the umbrella.
get_directory_property(hasParent PARENT_DIRECTORY)

# This adds the cmake directory to the module search paths,
# allowing us to use our modules.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# ----- Settings ------------------------------------------------------------ #
# This creates an ENABLE_TESTING option (default: ON) and enables the testing
# using the ctest executable.
# See: https://cmake.org/cmake/help/latest/manual/ctest.1.html
include(CTest)

# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

include(FindPackageHandleStandardArgs)

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_USE_WIN32_THREADS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
if (WIN32)
    # PThreads4WConfig.cmake given if `pthreads` installed via vcpkg..
    find_package(PThreads4W QUIET)
    # ... otherwise, overrides macros in Conda to find the `pthreads-win32`
    # package
    if (NOT PThreads4W_FOUND)

        find_path(PThreads4W_INCLUDE_DIR
                NAMES pthread.h
                PATHS $ENV{LIBRARY_INC}
                NO_DEFAULT_PATH
                DOC "PThreads4W include directory.")

        find_library(PThreads4W_LIBRARY
                NAMES pthreadVC2 pthreadVC3 pthread
                PATHS $ENV{LIBRARY_LIB}
                NO_DEFAULT_PATH
                DOC "PThreads4W library.")

        find_package_handle_standard_args(
                PThreads4W
                REQUIRED_VARS PThreads4W_LIBRARY PThreads4W_INCLUDE_DIR)
    endif ()

    if (PThreads4W_FOUND)
        if (NOT TARGET PThreads4W::PThreads4W)
            add_library(PThreads4W::PThreads4W UNKNOWN IMPORTED)
            set_target_properties(
                    PThreads4W::PThreads4W PROPERTIES
                    IMPORTED_LOCATION ${PThreads4W_LIBRARY}
                    INTERFACE_INCLUDE_DIRECTORIES ${PThreads4W_INCLUDE_DIR}
                    IMPORTED_LINK_INTERFACE_LANGUAGES "C")
        endif ()
    endif ()
endif ()

if (NOT hasParent)
    # If it's not being called by the umbrella, we need to
    # look for the system-installed libraries.
    find_package(MCFBlock REQUIRED)
endif ()

# ----- LEMON --------------------------------------------------------------- #
# abLEMON, the maintained fork of LEMON, is consumed as an installed package,
# found by cmake/FindLEMON.cmake as LEMON is, which defines the imported
# target LEMON::LEMON.
find_package(LEMON REQUIRED)

# ----- The warm start of abLEMON ------------------------------------------- #
# MCFLemonSolver re-optimizes with runWarm() of the four minimum cost flow
# algorithms, which abLEMON does not have: the build writes it into private
# copies of five headers of abLEMON, put on the include path before the
# installed ones [see shim/README.md]. abLEMON 1.9.0 is told from LEMON, and
# from the earlier abLEMON, by potentialMethod() in
# lemon/k_best_extreme_flows.h.
find_file(_ablemon_hdr lemon/k_best_extreme_flows.h
        PATHS ${LEMON_INCLUDE_DIRS} NO_DEFAULT_PATH)
if (_ablemon_hdr)
    file(STRINGS "${_ablemon_hdr}" _ablemon_pot REGEX "potentialMethod")
endif ()
if (NOT _ablemon_pot)
    message(FATAL_ERROR "MCFLemonSolver: the LEMON in ${LEMON_INCLUDE_DIRS} "
            "is not abLEMON 1.9.0 or later "
            "(https://github.com/dmeoli/abLEMON), which it needs")
endif ()
unset(_ablemon_hdr CACHE)
unset(_ablemon_pot)
set(LEMON_SHIM_DIR ${CMAKE_CURRENT_BINARY_DIR}/lemon-shim)
# the rewrites of shim/rewrites.txt, the same as shim/rewrite.awk of the
# makefiles: each anchor a whole line, found exactly once; the ';' of the
# rules are kept out of the way of the CMake lists
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/shim/rewrites.txt _shim_rules)
string(REPLACE ";" "<semicolon>" _shim_rules "${_shim_rules}")
string(REGEX REPLACE "#[^\n]*\n" "" _shim_rules "${_shim_rules}")
string(REPLACE "\n" ";" _shim_rules "${_shim_rules}")
list(REMOVE_ITEM _shim_rules "")
set(_shim_hdrs "")
foreach (_rule IN LISTS _shim_rules)
    string(FIND "${_rule}" "|" _p)
    string(SUBSTRING "${_rule}" 0 ${_p} _hdr)
    list(APPEND _shim_hdrs ${_hdr})
endforeach ()
list(REMOVE_DUPLICATES _shim_hdrs)
foreach (_hdr IN LISTS _shim_hdrs)
    unset(_hdr_src CACHE)
    find_file(_hdr_src lemon/${_hdr} PATHS ${LEMON_INCLUDE_DIRS}
            NO_DEFAULT_PATH)
    file(READ ${_hdr_src} _text)
    set(_text "\n${_text}")
    foreach (_rule IN LISTS _shim_rules)
        string(FIND "${_rule}" "|" _p)
        string(SUBSTRING "${_rule}" 0 ${_p} _rule_hdr)
        if (NOT _rule_hdr STREQUAL _hdr)
            continue ()
        endif ()
        math(EXPR _p "${_p} + 1")
        string(SUBSTRING "${_rule}" ${_p} -1 _rest)
        string(FIND "${_rest}" "|" _p)
        string(SUBSTRING "${_rest}" 0 ${_p} _anchor)
        math(EXPR _p "${_p} + 1")
        string(SUBSTRING "${_rest}" ${_p} 1 _op)
        math(EXPR _p "${_p} + 1")
        string(SUBSTRING "${_rest}" ${_p} -1 _arg)
        string(REPLACE "<semicolon>" ";" _anchor "${_anchor}")
        string(REPLACE "<semicolon>" ";" _arg "${_arg}")
        string(FIND "${_text}" "\n${_anchor}\n" _at)
        string(FIND "${_text}" "\n${_anchor}\n" _last REVERSE)
        if (_at EQUAL -1 OR NOT _at EQUAL _last)
            message(FATAL_ERROR "MCFLemonSolver: lemon/${_hdr} has not one "
                    "line \"${_anchor}\" [see shim/rewrites.txt]")
        endif ()
        if (_op STREQUAL "<")
            file(READ ${CMAKE_CURRENT_SOURCE_DIR}/shim/${_arg} _block)
            set(_new "\n${_block}${_anchor}\n")
        elseif (_op STREQUAL "+")
            set(_new "\n${_anchor}${_arg}\n")
        else ()
            set(_new "\n${_arg}\n")
        endif ()
        string(REPLACE "\n${_anchor}\n" "${_new}" _text "${_text}")
    endforeach ()
    string(SUBSTRING "${_text}" 1 -1 _text)
    file(WRITE ${LEMON_SHIM_DIR}/lemon/${_hdr} "${_text}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_hdr_src})
endforeach ()
unset(_hdr_src CACHE)
# a change of the rules or of the blocks writes the headers again
file(GLOB _shim_files ${CMAKE_CURRENT_SOURCE_DIR}/shim/*)
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_shim_files})

# Alias LEMON::LEMON to a plain `lemon` target name used throughout this file.
if (NOT TARGET lemon)
    add_library(lemon ALIAS LEMON::LEMON)
endif ()

# if SMSpp_which_insert_LEMON is *not* given in input as CMake option...
# (`NOT DEFINED SMSpp_which_insert_LEMON` tests `set( )`, not
# `add_compile_definitions( )`)
if (NOT DEFINED SMSpp_which_insert_LEMON)
    # ... then use 3 in MCFLemonSolver
    add_compile_definitions(SMSpp_which_insert_LEMON=3)
else ()
    add_compile_definitions(
            SMSpp_which_insert_LEMON=${SMSpp_which_insert_LEMON})
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- RPATH --------------------------------------------------------------- #
# The installed libraries and executables look for the libraries of the same
# installation relative to their own location, so that the installed tree can
# be moved as a whole; the directories of the external libraries are kept.
include(GNUInstallDirs)
if (NOT CMAKE_INSTALL_RPATH)
    if (APPLE)
        set(rpathOrigin "@loader_path")
    else ()
        set(rpathOrigin "$ORIGIN")
    endif ()
    file(RELATIVE_PATH binToLib ${CMAKE_INSTALL_FULL_BINDIR}
            ${CMAKE_INSTALL_FULL_LIBDIR})
    set(CMAKE_INSTALL_RPATH "${rpathOrigin}" "${rpathOrigin}/${binToLib}")
endif ()
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_20)

# The name of the shared library carries the version of its binary interface,
# which any minor release may change while the major version is 0, and then
# any major release.
if (PROJECT_VERSION_MAJOR EQUAL 0)
    set(soVersion ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
else ()
    set(soVersion ${PROJECT_VERSION_MAJOR})
endif ()
set_target_properties(${modName} PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${soVersion})

# The classes of the module register themselves in the factory from a static
# initialiser, and a linker that drops what looks unused takes the
# registration away with it: the target tells whoever links it to keep the
# symbol that forces the module in, or the library as a whole where naming
# the symbol is not enough.
if (BUILD_SHARED_LIBS)
    if (MSVC)
        target_link_options(${modName} INTERFACE
                "/INCLUDE:SMSpp_force_load_${modName}")
    elseif (UNIX AND (NOT APPLE))
        # the ELF linker decides by the symbols a regular object asks for, so
        # naming the symbol is not enough and the library has to be kept
        target_link_options(${modName} INTERFACE "LINKER:--no-as-needed")
    elseif (APPLE)
        # ld64 keeps a dylib that provides a symbol to the link, one asked
        # for by -u included, even under -dead_strip_dylibs
        target_link_options(${modName} INTERFACE
                "LINKER:-u,_SMSpp_force_load_${modName}")
    elseif (MINGW)
        # the PE linker imports a DLL for the symbols asked for, -u
        # included
        target_link_options(${modName} INTERFACE
                "LINKER:-u,SMSpp_force_load_${modName}")
    endif ()
endif ()

# Add the /bigobj flag for LEMON
if (MSVC)
    target_compile_options(${modName} PRIVATE /bigobj)
endif ()

# When adding source files with target_sources(), PRIVATE means that the files
# should only be added to this library, whereas PUBLIC means they should be
# added to this library and to any target that links to it.
# INTERFACE can be used for sources that should not be added to this library
# but should be added to anything that links to it.
# Note: do not GLOB files here.
target_sources(${modName} PRIVATE src/MCFLemonSolver.cpp)

# When using target_include_directories(), PUBLIC means that any targets
# that link to this target also need that include directory.
# Other options are PRIVATE (only affect the current target, not dependencies),
# and INTERFACE (only needed for dependencies).

# Different INSTALL_INTERFACE and BUILD_INTERFACE paths are used when
# generating the target import file (***Targets.cmake).
# This means that if a target finds this library in its build directory
# will look into the BUILD_INTERFACE path for its headers, if it finds it
# installed in the system will look into the INSTALL_INTERFACE path.
# LEMON's include directories are carried by the `lemon` target itself (see the
# LEMON section above), whether it comes from an installed package or from the
# fetched sources, so they no longer need to be repeated here.
target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

# The rewritten headers of abLEMON [see the warm start of abLEMON above] are
# found before the installed ones, by this library and by whatever is built
# against it, here or against the installed package, since MCFLemonSolver.h
# calls runWarm(); the installed copies are those of the abLEMON it was built
# with.
target_include_directories(${modName} BEFORE PUBLIC
        $<BUILD_INTERFACE:${LEMON_SHIM_DIR}>
        $<INSTALL_INTERFACE:include/${modNamespace}/lemon-shim>)

# When linking other targets to the library with target_link_libraries(),
# PRIVATE means that the libraries will be linked only to this library,
# PUBLIC means they will be linked also to the targets that depend on this
# library, INTERFACE means they will be linked only to the targets that depend
# on this library.
target_link_libraries(${modName} PUBLIC
        ${modNamespace}::MCFBlock
        lemon)

if (WIN32)
    target_link_libraries(${modName} PUBLIC PThreads4W::PThreads4W)
endif ()

# This alias is defined so that executables in this same project can use
# the library with this notation.
add_library(${modNamespace}::${modName} ALIAS ${modName})

# ----- C++20 module -------------------------------------------------------- #
# The module smspp.mcflemonsolver (src/MCFLemonSolver.cppm) exports the names
# of the headers and re-exports smspp.mcfblock; like the module smspp of the
# core, it is built when SMSpp_BUILD_MODULE is ON, which needs CMake 3.28 or
# later, the Ninja or Visual Studio generator and a compiler that scans the
# dependencies of the modules.
# Where these are missing, or a module it imports is not there, the option
# is turned off with a warning [see CheckModule.cmake of SMS++].
smspp_check_module(${modNamespace}::MCFBlock_module)
if (SMSpp_BUILD_MODULE)
    # the scan of the module dependencies, and $<LINK_ONLY> honoured also for
    # the target itself, which cmake_policy(VERSION 3.15) above leaves off
    cmake_policy(SET CMP0155 NEW)
    cmake_policy(SET CMP0131 NEW)

    add_library(${modName}_module)
    target_sources(${modName}_module PUBLIC
            FILE_SET CXX_MODULES
            BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/src
            FILES src/${modName}.cppm)
    target_compile_features(${modName}_module PUBLIC cxx_std_20)
    # the interface takes the include paths and the definitions of the
    # library, but not its compile options (e.g., the -fopenmp of a source,
    # which no header needs), so that it agrees with the module it imports
    target_include_directories(${modName}_module PRIVATE
            $<TARGET_PROPERTY:${modName},INTERFACE_INCLUDE_DIRECTORIES>)
    target_compile_definitions(${modName}_module PRIVATE
            $<TARGET_PROPERTY:${modName},INTERFACE_COMPILE_DEFINITIONS>)
    target_link_libraries(${modName}_module PUBLIC $<LINK_ONLY:${modName}>
            ${modNamespace}::MCFBlock_module)
    add_library(${modNamespace}::${modName}_module ALIAS ${modName}_module)
endif ()

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    add_subdirectory(test)
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # .dll
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} # .so .dylib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}) # .lib

# Install the module, whose interface is compiled again by each importer
if (SMSpp_BUILD_MODULE)
    install(TARGETS ${modName}_module
            EXPORT ${modName}Targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            FILE_SET CXX_MODULES
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}/modules)
endif ()

# Install the headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# Install the rewritten headers of abLEMON (see target_include_directories
# above).
install(DIRECTORY ${LEMON_SHIM_DIR}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace}/lemon-shim)

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# LEMON is consumed as an installed (imported) package, so it is neither
# installed nor exported here: it is already discoverable via
# find_package(LEMON) (see the find_dependency(LEMON) injected into the
# package config below).

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

# Downstream projects that do find_package(MCFLemonSolver) must also be able
# to locate LEMON: inject a find_dependency(LEMON) into the generated config,
# which is resolved there by the FindLEMON module we install alongside it
# (see below).
set(LEMON_FIND_DEPENDENCY "find_dependency(LEMON)")

# Generate the package configuration file, that allows other
# CMake projects to find the library with find_package().
configure_package_config_file(
        ${CMAKE_CURRENT_LIST_DIR}/cmake/${modName}Config.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the package version and configuration files.
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the custom find modules (e.g. FindLEMON) so that downstream projects
# can resolve the find_dependency() calls emitted in the package config above.
file(GLOB findMods ${CMAKE_CURRENT_LIST_DIR}/cmake/Find*)
install(FILES ${findMods}
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

# Install the README and LICENSE files.
install(FILES
        ${CMAKE_CURRENT_LIST_DIR}/README.md
        ${CMAKE_CURRENT_LIST_DIR}/LICENSE
        DESTINATION ${CMAKE_INSTALL_DATADIR}/${modName})

# ----- Add the build tree to RPATH------------------------------------------ #
# Add the binary directory to RPATH, so executables that use
# dynamic loading will look into it.
if (hasParent)
    set(SMSPP_BUILD_RPATH
            "${SMSPP_BUILD_RPATH};${CMAKE_CURRENT_BINARY_DIR}"
            PARENT_SCOPE)
endif ()

# ----- Packaging support --------------------------------------------------- #
# The following commands allow one to build a proper package using CPack.
# See: https://cmake.org/cmake/help/latest/manual/cpack.1.html

# Using this with the umbrella project is too messy.
if (NOT hasParent)
    set(CPACK_PACKAGE_VENDOR "SMS++ Team")
    set(CPACK_PACKAGE_CONTACT "Donato Meoli")
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY ${PROJECT_DESCRIPTION})
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_LIST_DIR}/LICENSE")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_LIST_DIR}/README.md")
    set(CPACK_SOURCE_IGNORE_FILES
            \\\\.git/
            \\\\.gitignore$
            \\\\.gitattributes$)
    include(CPack)
endif ()

# --------------------------------------------------------------------------- #
# Remove from the search paths the cmake directory we added at the beginning.
list(REMOVE_AT CMAKE_MODULE_PATH -1)

# --------------------------------------------------------------------------- #
