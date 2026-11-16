# --------------------------------------------------------------------------- #
#    CMake file for MCFLemonSolver/test                                       #
#                                                                             #
#    In this file we define the unit tests for this library.                  #
#    Tests that should be supported by ctest executable must be added         #
#    with the add_test() command.                                             #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Universita' di Pisa                             #
# --------------------------------------------------------------------------- #

# The classes of the modules register themselves in the factories from
# static objects, which nothing here names: a shared library is kept by
# what its module asks of the linker, while an archive has to be taken
# whole [see SMSppModule.cmake of SMS++].
smspp_keep_modules()

# ----- MCFLemonSolver_test ------------------------------------------------- #
add_executable(${modName}_test test.cpp)
target_link_libraries(${modName}_test PRIVATE SMS++::MCFLemonSolver)

smspp_copy_runtime(${modName}_test)

# The tester builds its instance in memory and needs no data; it writes its
# DMX file in the build directory and removes it.
add_test(NAME ${modName}_test
         COMMAND ${modName}_test
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ----- the label of this module's tests ------------------------------------ #

# The label lets the module CI select exactly the tests exercising this module
# with `ctest -L <module>` (see tests/cmake/TestLabels.cmake in the umbrella
# for the system-tests counterpart): without it these tests are built and
# never run.
get_property(_tests DIRECTORY PROPERTY TESTS)
if (_tests)
    set_tests_properties(${_tests} PROPERTIES LABELS "${modName}")
endif ()

# --------------------------------------------------------------------------- #
