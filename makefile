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
#   flows, are added to the factory                                          #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(SMS++INC)    = the -I$( core SMS++ directory )                 #
#           $(SMS++OBJ)    = the core SMS++ library                          #
#           $(MCFBkOBJ)    = the object(s) / library for MCFBlock            #
#           $(MCFBkINC)    = the -I$( source directory ) for MCFBlock        #
#           $(MCFLESDR)    = the directory where the source is               #
#                                                                            #
#   abLEMON is expected to be an installed package; its paths come from      #
#   extlib (LEMON_ROOT in extlib/makefile-default-paths-*). Its headers do   #
#   not have the warm start, which this makefile writes into private copies  #
#   of them that shadow the installed ones (see the README and the rule      #
#   below).                                                                  #
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

# LEMON (graph library) - - - - - - - - - - - - - - - - - - - - - - - - - - -
# LEMON paths come from extlib, exactly like the other external libraries:
# $(LEMON_ROOT) is set in extlib/makefile-default-paths-* (overridable via
# extlib/makefile-paths) and turned into $(libLEMONINC)/$(libLEMONLIB) here.
include $(MCFLESDR)/../extlib/makefile-libLEMON

# bare LEMON include dir (no -I, no quotes) used as the source for the shim
# rule; $(libLEMONBSCDIR) is the unquoted LEMON base dir exported by
# makefile-libLEMON
LEMON_INCDIR = $(libLEMONBSCDIR)/include

# The warm start: runWarm() of the four minimum cost flow algorithms, which
# abLEMON does not have, is written into private copies of five of its
# headers in $(MCFLESHIM), put on the include path before $(LEMON_INCDIR),
# by shim/rewrite.awk with the rules of shim/rewrites.txt, the same as the
# CMake build [see shim/README.md].
MCFLESHIM = $(MCFLESDR)/lemon-shim
LEMON_SHIM_HDRS = $(MCFLESHIM)/lemon/network_simplex.h \
                  $(MCFLESHIM)/lemon/cost_scaling.h \
                  $(MCFLESHIM)/lemon/cycle_canceling.h \
                  $(MCFLESHIM)/lemon/capacity_scaling.h \
                  $(MCFLESHIM)/lemon/bits/mcf_feasible_flow.h
LEMON_SHIM_INC = -I$(MCFLESHIM)

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFLEOBJ = $(MCFLESDR)/obj/MCFLemonSolver.o

MCFLEINC = -I$(MCFLESDR)/include $(LEMON_SHIM_INC) $(libLEMONINC)

MCFLEH   = $(MCFLESDR)/include/MCFLemonSolver.h

MCFLELIB = $(libLEMONLIB)

# the C++20 module smspp.mcflemonsolver - - - - - - - - - - - - - - - - - - -
# [see makefile-mod of SMS++]

$(eval $(call SMSppMODrule,MCFLE,$(MCFLESDR),MCFLemonSolver,\
	smspp.mcflemonsolver,MCFBk))

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

# abLEMON 1.9.0 is told from LEMON, and from the earlier abLEMON, by
# potentialMethod() in lemon/k_best_extreme_flows.h; a header is not written if
# an anchor of its rules is not there exactly once
LEMON_SHIM_SRC = $(MCFLESDR)/shim/rewrites.txt $(MCFLESDR)/shim/rewrite.awk \
	$(wildcard $(MCFLESDR)/shim/*.inc)

$(MCFLESHIM)/lemon/%.h: $(LEMON_INCDIR)/lemon/%.h $(LEMON_SHIM_SRC)
	@grep -qs potentialMethod \
	  $(LEMON_INCDIR)/lemon/k_best_extreme_flows.h || \
	  { echo "MCFLemonSolver: the LEMON in $(LEMON_INCDIR) is not" \
	    "abLEMON 1.9.0 or later, which it needs" >&2; exit 1; }
	mkdir -p $(dir $@)
	awk -v hdr=$*.h -v dir=$(MCFLESDR)/shim \
	  -f $(MCFLESDR)/shim/rewrite.awk $(MCFLESDR)/shim/rewrites.txt $< \
	  > $@ || { rm -f $@; exit 1; }

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MCFLESDR)/obj/MCFLemonSolver.o: $(MCFLESDR)/src/MCFLemonSolver.cpp \
	$(MCFLEH) $(LEMON_SHIM_HDRS) $(SMS++OBJ) $(MCFBkOBJ)
	$(CC) -c $(MCFLESDR)/src/MCFLemonSolver.cpp -o $@ $(SW) \
	$(SMS++INC) $(MCFBkINC) $(MCFLEINC) \
	-DSMSpp_which_insert_LEMON=$(which_insert_LEMON)

########################## End of makefile ###################################
