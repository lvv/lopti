#
#  MakeFile for newton_testset, subdirectory methods
#  Version 1.1
#  Latest Change     2nd June 2006
#  Written by: Lutz Weimann
#

ARCH = $(shell uname)

#CFLAGS = -g $(DFLAGS) -DDEBUG
CFLAGS = -O -DFMAT $(DFLAGS)
F1FLAGS = -O $(DFLAGS)

ifeq ($(ARCH),SunOS)
CC = cc
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

METHODS = nleq1.o nleq2.o linalg_nleq1.o linalg_nleq2.o \
          zibconst.o zibsec.o zibmon.o wnorm.o \
          nleq_err.o qnerr.o nleq_res.o  qnres.o utils.o \
          jacobian_and_linalg.o croutines_ftnint.o
          
$(LIBDIR)/libmethod.a: $(METHODS)
	ar r $(LIBDIR)/libmethod.a $(METHODS); \
	ranlib $(LIBDIR)/libmethod.a
	
clean: 
	rm -f $(METHODS)

cleanall:
	rm -f  $(METHODS) $(LIBDIR)/libmethod.a \
	croutines_ftnint.o nleq_err.o qnerr.o nleq_res.o qnres.o
	
utils.o jacobian_and_linalg.o: nleq.h jacobian.h
