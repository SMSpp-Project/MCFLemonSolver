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

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MCFLEOBJ = $(MCFLESDR)/obj/MCFLemonSolver.o

MCFLEINC = -I$(MCFLESDR)/include -I$(MCFLESDR)/lemon-development/lib/include

MCFLEH   = $(MCFLESDR)/include/MCFLemonSolver.h

MCFLELIB = -L$(MCFLESDR)/lemon-development/lib/lib -lemon

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFLEOBJ) $(MCFLESDR)/*~

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

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MCFLESDR)/obj/MCFLemonSolver.o: $(MCFLESDR)/src/MCFLemonSolver.cpp \
	$(MCFLEH) $(SMS++OBJ) $(MCFBkOBJ) 
	$(CC) -c $(MCFLESDR)/src/MCFLemonSolver.cpp -o $@ $(SW) \
	$(SMS++INC) $(MCFBkINC) $(MCFLEINC) \
	-DSMSpp_which_insert_LEMON=$(which_insert_LEMON)

########################## End of makefile ###################################
