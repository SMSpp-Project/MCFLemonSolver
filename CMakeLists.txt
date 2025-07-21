# --------------------------------------------------------------------------- #
#    Main CMake file for SMS++ LEMONSolver                                    #
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
cmake_minimum_required(VERSION 3.16)
cmake_policy(VERSION 3.15)
project(LEMONSolver
        VERSION 0.1.0
        DESCRIPTION "SMS++ LEMONSolver"
        HOMEPAGE_URL https://gitlab.com/smspp/lemonsolver
        LANGUAGES C CXX)

# These variables make the code harder to read but easier to change.
set(modName ${PROJECT_NAME})
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
#include(CTest)

# Sets the default build type (if none was specified).
# See: https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html
include(BuildType)

# ----- Requirements -------------------------------------------------------- #
# This sets the variable CMAKE_USE_WIN32_THREADS_INIT, see:
# https://cmake.org/cmake/help/latest/module/FindThreads.html
if (WIN32)
    # PThreads4WConfig.cmake given if `pthreads` installed via vcpkg..
    find_package(PThreads4W QUIET)
    # ... otherwise, overrides macros in Conda to find the `pthreads-win32` package
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

        if (PThreads4W_INCLUDE_DIR AND PThreads4W_LIBRARY)
            add_library(PThreads4W::PThreads4W UNKNOWN IMPORTED)
            set_target_properties(
                    PThreads4W::PThreads4W PROPERTIES
                    IMPORTED_LOCATION ${PThreads4W_LIBRAY}
                    IMPORTED_LOCATION_DEBUG ${PThreads4W_LIBRAY}
                    INTERFACE_INCLUDE_DIRECTORIES ${PThreads4W_INCLUDE_DIR})
            set(PThreads4W_FOUND TRUE)
        else ()
            set(_missing "")
            if (NOT PThreads4W_INCLUDE_DIR)
                string(APPEND _missing " PThreads4W_INCLUDE_DIR")
            endif ()
            if (NOT PThreads4W_LIBRARY)
                string(APPEND _missing " PThreads4W_LIBRARY")
            endif ()
            message(FATAL_ERROR "Could NOT find PThreads4W (missing:${_missing})")
        endif ()
    endif ()
endif ()

if (NOT hasParent)
    # If it's not being called by the umbrella, we need to
    # look for the system-installed libraries.
    find_package(MCFBlock REQUIRED)
endif ()

# Define a function to wrap add_subdirectory with CMAKE_SUPPRESS_DEVELOPER_WARNINGS, i.e., -Wno-dev, set to TRUE
function(add_subdirectory_no_dev_warnings dir)
    set(old_value ${CMAKE_SUPPRESS_DEVELOPER_WARNINGS})
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS TRUE CACHE INTERNAL "Suppress developer warnings")
    add_subdirectory(${dir})
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ${old_value} CACHE INTERNAL "Suppress developer warnings")
endfunction()

# Set the lemon-development options in the cache
set(LEMON_ENABLE_GLPK OFF CACHE BOOL "Enable GLPK solver backend.")
set(LEMON_ENABLE_ILOG OFF CACHE BOOL "Enable ILOG (CPLEX) solver backend.")
set(LEMON_ENABLE_COIN OFF CACHE BOOL "Enable COIN solver backend.")
set(LEMON_ENABLE_SOPLEX OFF CACHE BOOL "Enable SoPlex solver backend.")

# Use the function to add the subdirectory
add_subdirectory_no_dev_warnings(lemon-development)

# if SMSpp_which_insert_LEMON is *not* given in input as CMake option...
# (`NOT DEFINED SMSpp_which_insert_LEMON` tests `set( )`, not `add_compile_definitions( )`)
if (NOT DEFINED SMSpp_which_insert_LEMON)
    # ... then use 3 in LEMONSolver
    add_compile_definitions(SMSpp_which_insert_LEMON=3)
else ()
    add_compile_definitions(SMSpp_which_insert_LEMON=${SMSpp_which_insert_LEMON})
endif ()

# ----- Configuration header ------------------------------------------------ #
# This will generate a *Config.h header in the build directory.
configure_file(cmake/${modName}Config.h.in ${modName}Config.h)

# ----- Library ------------------------------------------------------------- #
# With the BUILD_SHARED_LIBS variable we can specify if the library will
# be STATIC or SHARED, so no reason to do it now.
add_library(${modName})
target_compile_features(${modName} PUBLIC cxx_std_20)

# Add the /bigobj flag for lemon-development
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
target_include_directories(
        ${modName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/${modNamespace}>)

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

# ----- Subdirectories ------------------------------------------------------ #
if (BUILD_TESTING)
    # Just consuming BUILD_TESTING to avoid warnings
endif ()

# ----- Install instructions ------------------------------------------------ #
# The following commands are used when installing the library
# and its CMake configuration files on the system.
# They are not required for local builds (see below).

include(GNUInstallDirs)

# Install the library
install(TARGETS ${modName}
        EXPORT ${modName}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install the headers
install(DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# Install the auto-generated configuration header (see above).
install(FILES ${PROJECT_BINARY_DIR}/${modName}Config.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${modNamespace})

# Install the lemon library as well
install(TARGETS lemon
        EXPORT lemonTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# (Generate and) install the target import file, that allows other
# CMake projects to import the target.
install(EXPORT ${modName}Targets
        NAMESPACE ${modNamespace}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${modName})

install(EXPORT lemonTargets
        NAMESPACE lemon::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lemon)

# Generate the package version file, that allows other
# CMake projects to know the version.
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/${modName}ConfigVersion.cmake
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY AnyNewerVersion)

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
