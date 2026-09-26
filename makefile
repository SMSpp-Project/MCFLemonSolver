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
#   LEMON is expected to be an installed package; its paths come from        #
#   extlib (LEMON_ROOT in extlib/makefile-default-paths-*). The packaged     #
#   LEMON needs C++20 and MSVC portability fixes, so this makefile rewrites  #
#   the affected headers into a private shim that shadows them (see the      #
#   README and the rule below).                                              #
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

# bare LEMON include dir (no -I, no quotes) used as the source for the shim rule;
# $(libLEMONBSCDIR) is the unquoted LEMON base dir exported by makefile-libLEMON
LEMON_INCDIR = $(libLEMONBSCDIR)/include

# Compatibility shim: the packaged LEMON needs a C++20 fix (array_map.h, path.h
# still call the removed std::allocator::construct/destroy), an MSVC fix
# (adaptors.h LEMON_SCOPE_FIX), a missing-include fix (capacity_scaling.h uses
# RangeMap but some packagings, e.g. MacPorts, omit #include <lemon/maps.h>) and
# a numerical one (network_simplex.h takes the sign of a reduced cost, and the
# flow left on its artificial arcs, exactly, which with fractional costs makes
# it pivot for ever or report infeasibility). We
# rewrite private copies of those headers (see the rules below) and put
# $(MCFLESHIM) on the include path *before* $(LEMON_INCDIR) so they shadow the
# system ones. Mirrors the CMake build; all rewrites are idempotent and harmless
# on non-matching headers.
MCFLESHIM = $(MCFLESDR)/lemon-cxx20-shim
LEMON_SHIM_HDRS = $(MCFLESHIM)/lemon/bits/array_map.h \
                  $(MCFLESHIM)/lemon/path.h \
                  $(MCFLESHIM)/lemon/adaptors.h \
                  $(MCFLESHIM)/lemon/capacity_scaling.h \
                  $(MCFLESHIM)/lemon/network_simplex.h

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFLEOBJ = $(MCFLESDR)/obj/MCFLemonSolver.o

MCFLEINC = -I$(MCFLESDR)/include -I$(MCFLESHIM) $(libLEMONINC)

MCFLEH   = $(MCFLESDR)/include/MCFLemonSolver.h

MCFLELIB = $(libLEMONLIB)

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
# NOTE: '#' is escaped as '\#' because this is a make *variable* (unlike a recipe
# line, an unescaped '#' would start a make comment and truncate the value).
LEMON_SHIM_SED = \
	  -e 's/([A-Za-z_][A-Za-z0-9_]*)\.construct\(/std::allocator_traits<decltype(\1)>::construct(\1, /g' \
	  -e 's/([A-Za-z_][A-Za-z0-9_]*)\.destroy\(/std::allocator_traits<decltype(\1)>::destroy(\1, /g' \
	  -e 's/\#define LEMON_SCOPE_FIX\(OUTER, NESTED\) OUTER::NESTED/\#define LEMON_SCOPE_FIX(OUTER, NESTED) typename OUTER::template NESTED/'

$(MCFLESHIM)/lemon/%.h: $(LEMON_INCDIR)/lemon/%.h
	mkdir -p $(dir $@)
	sed -E $(LEMON_SHIM_SED) $< > $@

# capacity_scaling.h needs the base rewrites plus the maps.h include (RangeMap);
# this explicit rule overrides the generic one above for that single header. The
# extra expression adds #include <lemon/maps.h> right after #include <lemon/core.h>.
$(MCFLESHIM)/lemon/capacity_scaling.h: $(LEMON_INCDIR)/lemon/capacity_scaling.h
	mkdir -p $(dir $@)
	sed -E $(LEMON_SHIM_SED) \
	  -e '/#include <lemon\/maps.h>/d' \
	  -e 's@(#include <lemon/core.h>)@\1\n#include <lemon/maps.h>@' \
	  $< > $@

# network_simplex.h needs the base rewrites plus the two of its own, the same
# as the CMake build: the tolerances of the pivot rules and of the feasibility
# check, zero for an exact type so that on integer data the algorithm is
# unchanged, and the warm start of the basis. The code of both lives in shim/
# and sed inserts it at its anchor, so that it is written once and read as
# C++ [see shim/README.md].

$(MCFLESHIM)/lemon/network_simplex.h: $(LEMON_INCDIR)/lemon/network_simplex.h \
	  $(MCFLESDR)/shim/ns_eps.inc $(MCFLESDR)/shim/ns_runwarm.inc
	mkdir -p $(dir $@)
	sed -E $(LEMON_SHIM_SED) \
	  -e 's/^      int _search_arc_num;$$/      int _search_arc_num; Cost _eps;/' \
	  -e 's/_search_arc_num\(ns\._search_arc_num\)/&, _eps(ns_pivot_eps(ns._cost, ns._arc_num))/' \
	  -e 's/if \(c < 0\)/if (c < -_eps)/g' \
	  -e 's/min < 0/min < -_eps/g' \
	  -e 's/if \(min >= 0\) return/if (min >= -_eps) return/' \
	  -e 's/else if \(c >= 0\) \{/else if (c >= -_eps) {/' \
	  -e 's@^      // Check feasibility@      const Value feps = ns_flow_eps(_supply, _node_num); // Check feasibility@' \
	  -e 's/if \(_flow\[e\] != 0\) return/if (_flow[e] > feps || _flow[e] < -feps) return/' \
	  -e 's/^    int _root;$$/    int _root; bool _warm; bool _repair; int _unb_arc;/' \
	  -e 's@^      // Check the number types@      _warm = false; _repair = false; _unb_arc = -1; // Check the number types@' \
	  -e 's/if \(!initialPivots\(\)\) return UNBOUNDED;/if (!initialPivots()) { _unb_arc = in_arc; return UNBOUNDED; }/' \
	  -e 's/if \(delta >= MAX\) return UNBOUNDED;/if (delta >= MAX) { _unb_arc = in_arc; return UNBOUNDED; }/' \
	  -e 's@^        _upper\[_arc_id\[a\]\] = map\[a\];$$@        warmUpper(_arc_id[a], map[a]);@' \
	  -e 's@^        _supply\[_node_id\[n\]\] = map\[n\];$$@        warmSupply(_node_id[n], map[n]);@' \
	  -e 's/^(    NetworkSimplex\& lowerMap\(const .*\& map\) \{)$$/\1 _warm = false;/' \
	  -e 's/^(    NetworkSimplex\& stSupply\(const Node\& s, const Node\& t, Value k\) \{)$$/\1 _warm = false;/' \
	  -e 's/^(    NetworkSimplex\& reset\(\) \{)$$/\1 _warm = false;/' \
	  -e '/^namespace lemon \{$$/r $(MCFLESDR)/shim/ns_eps.inc' \
	  -e '/^      return start\(pivot_rule\);$$/{r $(MCFLESDR)/shim/ns_runwarm.inc' \
	  -e 'd}' \
	  $< > $@

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MCFLESDR)/obj/MCFLemonSolver.o: $(MCFLESDR)/src/MCFLemonSolver.cpp \
	$(MCFLEH) $(LEMON_SHIM_HDRS) $(SMS++OBJ) $(MCFBkOBJ)
	$(CC) -c $(MCFLESDR)/src/MCFLemonSolver.cpp -o $@ $(SW) \
	$(SMS++INC) $(MCFBkINC) $(MCFLEINC) \
	-DSMSpp_which_insert_LEMON=$(which_insert_LEMON)

########################## End of makefile ###################################
