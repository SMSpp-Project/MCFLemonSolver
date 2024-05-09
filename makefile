##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MCFLEMONSolver                                               #
#                                                                            #
#   The makefile takes in input the -I directives for all the external       #
#   libraries needed by MCFLEMONSolver, i.e., core SMS++ and MCFBlock.       #
#   These are *not* copied into $(MCFLEINC): adding those -I directives to   #
#   the compile commands will have to done by whatever "main" makefile is    #
#   using this. Analogously, any external library and the corresponding      #
#   -L< libdirs > will have to be added to the final linking command by      #
#   whatever "main" makefile is using this.                                  #
#                                                                            #
#   Note that, conversely, $(SMS++INC) is also assumed to include any        #
#   -I directive corresponding to external libraries needed by SMS++, at     #
#   least to the extent in which they are needed by the parts of SMS++       #
#   used by MCFLEMONSolver.                                                  #
#                                                                            #
#   Input:  $(CC)          = compiler command                                #
#           $(SW)          = compiler options                                #
#           $(MCFLESDR)    = the directory where the source is               #
#                                                                            #
#   Output: $(MCFLEOBJ)    = the final object(s) / library                   #
#           $(MCFLEH)      = the .h files to include                         #
#           $(MCFLEINC)    = the -I$( source directory )                     #
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

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MCFLEOBJ) $(MCFLESDR)/*~

# dependencies: every .o from its .cpp + every recursively included .h- - - -

$(MCFLESDR)/obj/MCFLemonSolver.o: $(MCFLESDR)/src/MCFLemonSolver.cpp \
	$(MCFLEH)
	$(CC) -c $(MCFLESDR)/src/MCFLemonSolver.cpp -o $@ $(SW) \
	$(MCFLEINC) 

########################## End of makefile ###################################
