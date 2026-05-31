##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MCFLemonSolver                                               #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by MCFLemonSolver, i.e., core SMS++ and MCFBlock.       #
#   These are *not* copied into $(MCFLEINC): adding those -I directives to   #
#   the compile commands will have to done by whatever "main" makefile is    #
#   using this. Analogously, any external library and the corresponding      #
#   -L< libdirs > will have to be added to the final linking command by      #
#   whatever "main" makefile is using this.                                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by MCFLemonSolver.                                                  #
#                                                                            #
#   This makefile internally defines the crucial macro which_insert_LEMON    #
#   which is passed in the compilation of MCFLemonSolver.cpp as              #
#   SMSpp_which_insert_LEMON and whose value, coded bit-wise, controls       #
#   which variants of the algorithm, in terms of allowed types of costs and  #
#   flows, are added to the factory
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(MCFBkOBJ)    = the object(s) / library for MCFBlock            #
#           $(MCFBkINC)    = the -I$( source directory ) for MCFBlock        #
#           $(MCFLESDR)    = the directory where the source is               #
#                                                                            #
#   LEMON is expected to be an installed package (no in-tree clone);         #
#   override LEMON_ROOT below if it lives elsewhere. Since the packaged      #
#   LEMON predates C++20, this makefile rewrites the two headers using       #
#   the removed std::allocator::construct/destroy into a private shim        #
#   that shadows them (see the README and the rule below).                   #
#                                                                            #
#   Output: $(MCFLEOBJ)    = the final object(s) / library                   #
#           $(MCFLEH)      = the .h files to include                         #
#           $(MCFLEINC)    = the -I$( source directory )                     #
#           $(MCFLELIB)    = external libraries + -L< libdirs >              #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# LEMON (graph library) paths - - - - - - - - - - - - - - - - - - - - - - - -
# LEMON is an installed package; override LEMON_ROOT (or LEMON_INCDIR /
# LEMON_LIBDIR) if it lives outside the default /usr prefix.
LEMON_ROOT   ?= /usr
LEMON_INCDIR ?= $(LEMON_ROOT)/include
LEMON_LIBDIR ?= $(LEMON_ROOT)/lib

# Compatibility shim: the packaged LEMON needs a C++20 fix (array_map.h, path.h
# still call the removed std::allocator::construct/destroy) and an MSVC fix
# (adaptors.h LEMON_SCOPE_FIX). We rewrite private copies of those headers (see
# the rule below) and put $(MCFLESHIM) on the include path *before* $(LEMON_INCDIR)
# so they shadow the system ones. Mirrors the CMake build; both rewrites are
# idempotent and harmless on non-matching headers.
MCFLESHIM = $(MCFLESDR)/lemon-cxx20-shim
LEMON_SHIM_HDRS = $(MCFLESHIM)/lemon/bits/array_map.h \
                  $(MCFLESHIM)/lemon/path.h \
                  $(MCFLESHIM)/lemon/adaptors.h

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFLEOBJ = $(MCFLESDR)/obj/MCFLemonSolver.o

MCFLEINC = -I$(MCFLESDR)/include -I$(MCFLESHIM) -I$(LEMON_INCDIR)

MCFLEH   = $(MCFLESDR)/include/MCFLemonSolver.h

MCFLELIB = -L$(LEMON_LIBDIR) -llemon

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFLEOBJ) $(MCFLESDR)/*~
	rm -rf $(MCFLESHIM)

# internal macros - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# the value of this macro is passed as SMSpp_which_insert_LEMON when
# MCFLemonSolver.cpp is compiled and controls which types are allowed for
# flows and costs:
# - bit 0 (+1): double
# - bit 1 (+2): long
# - bit 2 (+4): int
# note that 1 bit = 8 variants, 2 bits = 32 variants, 3 bits = 96 variants,
# hence the compile time may increase considerably with more bits on

which_insert_LEMON = 3

# Compatibility shim: rewrite the offending LEMON headers into $(MCFLESHIM) —
# the C++20 "<alloc>.construct(...)"/"<alloc>.destroy(...)" calls into the
# std::allocator_traits<decltype(<alloc>)>::construct/destroy form, and the MSVC
# LEMON_SCOPE_FIX definition into its portable form. Each rewrite is idempotent
# and harmless on non-matching headers, matching the CMake build.
$(MCFLESHIM)/lemon/%.h: $(LEMON_INCDIR)/lemon/%.h
	mkdir -p $(dir $@)
	sed -E \
	  -e 's/([A-Za-z_][A-Za-z0-9_]*)\.construct\(/std::allocator_traits<decltype(\1)>::construct(\1, /g' \
	  -e 's/([A-Za-z_][A-Za-z0-9_]*)\.destroy\(/std::allocator_traits<decltype(\1)>::destroy(\1, /g' \
	  -e 's/#define LEMON_SCOPE_FIX\(OUTER, NESTED\) OUTER::NESTED/#define LEMON_SCOPE_FIX(OUTER, NESTED) typename OUTER::template NESTED/' \
	  $< > $@

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MCFLESDR)/obj/MCFLemonSolver.o: $(MCFLESDR)/src/MCFLemonSolver.cpp \
	$(MCFLEH) $(LEMON_SHIM_HDRS) $(SMS++OBJ) $(MCFBkOBJ)
	$(CC) -c $(MCFLESDR)/src/MCFLemonSolver.cpp -o $@ $(SW) \
	$(SMS++INC) $(MCFBkINC) $(MCFLEINC) \
	-DSMSpp_which_insert_LEMON=$(which_insert_LEMON)

########################## End of makefile ###################################
