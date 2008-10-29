#
#  MakeFile for newton_testset, subdirectory problems
#  Version 1.1
#  Latest Change     2nd June 2006
#  Written by: Lutz Weimann
#

ARCH = $(shell uname)

F1FLAGS = -O $(DFLAGS)

ifeq ($(ARCH),SunOS)
#LDFLAGS = -Bstatic
FFLAGS = $(F1FLAGS)
endif    #ifeq ($(ARCH),SunOS)

ifeq ($(ARCH),Darwin)
ifeq ($(FC),f77)
FFLAGS = -qextname $(F1FLAGS)
endif    #ifeq ($(FC),f77)
ifeq ($(FC),g77)
FFLAGS = -fno-second-underscore $(F1FLAGS)
endif    #ifeq ($(FC),g77)
endif    #ifeq ($(ARCH),Darwin)

ifeq ($(ARCH),Linux)
FFLAGS = -fno-second-underscore $(F1FLAGS)
endif    #ifeq ($(ARCH),Linux)

LIBDIR = ../lib

PROBLEMS = CompFacRK.o \
          CREquiConv4.o \
          CRSteadyState5.o \
          Cutlip.o \
          dcpffj.o \
          dcprfj.o \
          delhfj.o \
          destfj.o \
          dfdcfj.o \
          dficfj.o \
          dgupfj.o \
          dhhdfj.o \
          dierfj.o \
          disofj.o \
          dmetfj.o \
          dsfdfj.o \
          dsfifj.o \
          dsulfj.o \
          BrownAlmLin.o \
          BroydenBanded.o \
          BroydenTridiag.o \
          Chebyquad.o \
          DiscreteBoundary.o \
          DiscreteIntegral.o \
          HelicalValley.o \
          PowellBadlyScaled.o \
          PowellSingular.o \
          Rosenbrock.o \
          Trigonometric.o \
          VariablyDim.o \
          Watson.o \
          Wood.o \
          ChemEqui2.o \
          LargeExp.o \
          ExpSin.o \
          ChemEqui3.o \
          ChemEqui1.o \
          DistillColumn.o \
          sst.o \
          HodgkinHuxley1.o \
          InfReflux.o \
          Order10to11.o \
          PartMethaneOx.o \
          ReinforcedBeam.o \
          Wallis1685.o \
          EsterificRea.o \
          Gupta.o \
          DewPoint.o \
          PartOxiOfMet.o \
          SulpDiox.o \
          InterSEllHy.o \
          faces.o
          
$(LIBDIR)/libproblem.a: $(PROBLEMS)
	ar r $(LIBDIR)/libproblem.a $(PROBLEMS); \
	ranlib $(LIBDIR)/libproblem.a
	
clean: 
	rm -f $(PROBLEMS)

cleanall:
	rm -f  $(PROBLEMS) $(LIBDIR)/libproblem.a
