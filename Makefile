 # $Id: Makefile 227 2010-11-25 14:30:08Z krasznaa $
 ###########################################################################
 # @Project: SFrame - ROOT-based analysis framework for ATLAS              #
 # @Package: Core                                                          #
 #                                                                         #
 # @author Stefan Ask       <Stefan.Ask@cern.ch>           - Manchester    #
 # @author David Berge      <David.Berge@cern.ch>          - CERN          #
 # @author Johannes Haller  <Johannes.Haller@cern.ch>      - Hamburg       #
 # @author A. Krasznahorkay <Attila.Krasznahorkay@cern.ch> - CERN/Debrecen #
 #                                                                         #
 # Makefile compiling the SFrameCore library.                              #
 #                                                                         #
 ###########################################################################

# Makefile updated by Jared Vasquez on Jan 29th, 2017
# Package information
APP = quickCLs
LIBRARY = quickCLs
OBJDIR  = obj
DEPDIR  = $(OBJDIR)/dep
SRCDIR  = src
INCDIR  = inc

ifeq ($(shell root-config --platform), macosx)
	BOOSTLIBS = -lboost_system-mt -lboost_filesystem-mt -lboost_program_options-mt -lboost_regex-mt
else
# Obsolete (HSG7 ROOT 5)
	# BOOSTLIBS = -L/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/boost/boost-1.55.0-python2.7-x86_64-slc6-gcc48/boost-1.55.0-python2.7-x86_64-slc6-gcc48/lib/ -lboost_filesystem -lboost_program_options -lboost_regex -lboost_system
	# BOOST_INC = -I/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/boost/boost-1.55.0-python2.7-x86_64-slc6-gcc48/boost-1.55.0-python2.7-x86_64-slc6-gcc48/include/

	BOOSTLIBS = -L/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/boost/boost-1.60.0-python2.7-x86_64-slc6-gcc49/boost-1.60.0-python2.7-x86_64-slc6-gcc49/lib/ -lboost_filesystem -lboost_program_options -lboost_regex -lboost_system
	BOOST_INC = -I/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/boost/boost-1.60.0-python2.7-x86_64-slc6-gcc49/boost-1.60.0-python2.7-x86_64-slc6-gcc49/include/
	INCLUDES += $(BOOST_INC) 
endif


MOREROOTLIBS = -lRooFit -lRooFitCore -lRooStats -lHistFactory -lMinuit -lMathMore -lSmatrix


# Overwrite the default rule defined in Makefile.common
coredefault: default $(_BIN_PATH)/$(APP)

# Include the library compilation rules
include $(_DIRCLS)/Makefile.common

#
# Rules for compiling the executable
#
# Reminder: $(ROOTLIBS) comes from Makefile.arch which is picked up from the ROOT
# sources by Makefile.common...
#

$(_BIN_PATH)/$(APP): $(APP).o $(SHLIBFILE)
	@echo "Linking " $@
	@$(LD) $(LDFLAGS) -O2 $(OBJDIR)/$(APP).o -L$(_LIB_PATH) -l$(LIBRARY) \
		$(ROOTLIBS) $(MOREROOTLIBS) $(BOOSTLIBS) -o $@

$(APP).o: app/$(APP).cxx 
	@echo "Compiling $<"
	@mkdir -p $(OBJDIR)
	@$(CXX) $(CXXFLAGS) -O2 -c $< -o $(OBJDIR)/$(notdir $@) $(INCLUDES)
