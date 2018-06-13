
#To compile with gsl, uncomment the following line
USE_GSL=-DGSL

ifdef USE_GSL
GSLLIB=-lgsl
else
GSLLIB=
endif

CC=mpic++
CPPFLAGS= -w -O3 -c $(USE_GSL) $(INCLUDES)
SRCS=   TaxonSet.cpp Tree.cpp Random.cpp SequenceAlignment.cpp CodonSequenceAlignment.cpp \
	StateSpace.cpp CodonStateSpace.cpp ZippedSequenceAlignment.cpp SubMatrix.cpp \
	linalg.cpp Chrono.cpp \
	ProfileProcess.cpp \
	OneProfileProcess.cpp \
	SubstitutionProcess.cpp \
	PoissonSubstitutionProcess.cpp \
	PhyloProcess.cpp PoissonPhyloProcess.cpp \
	Bipartition.cpp BipartitionList.cpp TaxaParameters.cpp TreeList.cpp PolyNode.cpp correl.cpp correlation.cpp NNI.cpp\
        GammaBranchProcessVI.cpp GammaRateProcessVI.cpp  PoissonDPProfileProcessVI.cpp PoissonMixtureProfileProcessVI.cpp PoissonSBDPProfileProcessVI.cpp RASCATGammaPhyloProcessVI.cpp RASCATSBDPGammaPhyloProcessVI.cpp \
        MixtureProfileProcessVI.cpp DPProfileProcessVI.cpp SBDPProfileProcessVI.cpp RateProcessVI.cpp BranchProcessVI.cpp\


OBJS=$(patsubst %.cpp,%.o,$(patsubst %.c,%.o,$(SRCS)))
ALL_SRCS:=$(wildcard *.cpp *.c)
ALL_OBJS=$(patsubst %.cpp,%.o,$(patsubst %.c,%.o,$(ALL_SRCS)))

PROGSDIR=../data
ALL= pbVI_mpi
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$


# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

$(PROGSDIR)/pbVI_mpi: PBVI.o $(OBJS)
	$(CC) PBVI.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d
	-rm -f $(PROGS)

