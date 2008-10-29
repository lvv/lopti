#
#  MakeFile for newton_testset
#  Version 1.1
#  Latest Change     2nd June 2006
#  Written by: Lutz Weimann
#

ARCH = $(shell uname)

#CFLAGS = -g $(DFLAGS) -DDEBUG
CFLAGS = -O $(DFLAGS)
F1FLAGS = -O $(DFLAGS)

ifeq ($(ARCH),SunOS)
CC = cc
#VECLIB = -lsunperf -lsunmath -lm -lF77 -lM77
#VECLIB = -lsunperf -lsunmath -lm
LD = $(FC)
#LDFLAGS = -Bstatic
FFLAGS = $(F1FLAGS)
endif    #ifeq ($(ARCH),SunOS)

ifeq ($(ARCH),Darwin)
#VECLIB = -framework vecLib
LD = $(FC)
ifeq ($(FC),f77)
FFLAGS = -qextname $(F1FLAGS)
ADDLIB = -lSystemStubs
endif    #ifeq ($(FC),f77)
ifeq ($(FC),g77)
FFLAGS = -fno-second-underscore $(F1FLAGS)
ADDLIB = -lSystemStubs -L/usr/lib/gcc/powerpc-apple-darwin8/4.0.1 -lgcc
endif    #ifeq ($(FC),g77)
endif    #ifeq ($(ARCH),Darwin)

ifeq ($(ARCH),Linux)
#VECLIB = -llapack -lblas
LD = $(FC)
FFLAGS = -fno-second-underscore $(F1FLAGS)
endif    #ifeq ($(ARCH),Linux)

VECLIB = lapack_and_blas.o

PRODIR := $(PWD)/problems
METHDIR := $(PWD)/methods
LIBDIR := $(PWD)/lib

TESTSET = testset.o

testset:  $(TESTSET) $(LIBDIR) $(LIBDIR)/libproblem.a $(LIBDIR)/libmethod.a \
          $(VECLIB)
	$(LD) $(TESTSET) -o testset -L$(LIBDIR) -lproblem -lmethod \
	      $(VECLIB) $(ADDLIB)

$(LIBDIR): 
	mkdir $(LIBDIR)

$(LIBDIR)/libproblem.a: now
	cd $(PRODIR); make

$(LIBDIR)/libmethod.a: now
	cd $(METHDIR); make

clean:
	cd $(PRODIR) ; make clean
	cd $(METHDIR) ; make clean
	rm -f $(TESTSET) lapack_and_blas.o

cleanall:
	cd $(PRODIR) ; make cleanall
	cd $(METHDIR) ; make cleanall
	rm -f $(TESTSET) lapack_and_blas.o testset

now:
