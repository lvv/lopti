
# this include is part of lvvlib
include  ../lvv/include.mk
	
SPEED ?= DEBUG
PREFIX ?= /usr/local
VERSION = 0.3
#FCFLAGS +=  -O2
FCFLAGS += -frange-check -fbounds-check -O0 -ggdb3
LDFLAGS += -L. -L /usr/local/lib  -lgsl -lgslcblas -lcondor -llopti -lgfortran 
CXXFLAGS += -std=c++0x -I /usr/local/include/ -I .. -I .
XGRAPHIC = xgraphic  -scat -markcol=-1  -g2 -logy -leg -legpos=3  -legsiz=1 -legtyp=2 -titgen="Convergance speed for dirivative-free algorithms" -titx="Objective function evaluation count" -tity="Distance to optimum:  log10 ( | X - X_opt | )" 

#.DEFAULT_GOAL := t-lopti-r
#.PHONY: t-lopti-r all

test:	t-lopti
	./t-lopti

install: liblopti.so
	mkdir -p		$(PREFIX)/include/lopti/
	cp -v *.h 		$(PREFIX)/include/lopti/
	mkdir -p		$(PREFIX)/lib/
	cp -va liblopti.so*	$(PREFIX)/lib/

clean:
	rm -f liblopti.so*
	rm -f {,newuoa/}*.o
	rm -f doc/*.{html,css}
	rm -f t-{lopti,newuoa}
	rm -f gmon.out

 
liblopti.so: *.h newuoa/bigden.o newuoa/biglag.o newuoa/calfun.o newuoa/trsapp.o newuoa/update.o 
	$(FC) $(FCFLAGS) -fPIC -c newuoa/{bigden,biglag,calfun,trsapp,update}.f
	gcc -shared -Wl,-soname,liblopti.so.$(VERSION) -o liblopti.so.$(VERSION) *.o
	ln -sf liblopti.so.$(VERSION) liblopti.so

t-%: t-%.cc *.h liblopti.so
t-%: SPEED=DEBUG

t-lopti: t-lopti.cc liblopti.so

t-lopti-xg: t-lopti log/condor
	$<
	cd log; LANG=C.iso88591 $(XGRAPHIC) `find  -type f \( -newer ../$< -o -cmin -1 \) -printf "%f\n"`

log/condor: t-lopti 

t-lopti-ps: t-lopti.ps
t-lopti.ps: t-lopti t-lopti-r 
	$<
	cd log; $(XGRAPHIC) -pscoltex=../$@  `find  -type f \( -newer ../$< -o -cmin -1 \) -printf "%f\n"`
	gsview $@

t-condor: CXXFLAGS +=   -I ..
t-condor: LDFLAGS  +=   -L /usr/local/lib/ -lcondor  -lm 

t-condor: t-condor.cc condor-wrap.h
	$(CXX) $(CXXFLAGS) -DOPTI=CONDOR  $< -o $@ $(LDFLAGS)

t-nm: t.cc  gsl-nelder-mead-wrap.h
	$(CXX) $(CXXFLAGS) -DOPTI=NM  $< -o $@ $(LDFLAGS)

t-newuoa-old: newuoa/t-newuoa.cc newuoa-wrap.h  lopti.h
	#	rm -f *.o
	#$(FC)  $(FCFLAGS) -c bigden.f  biglag.f  calfun.f   trsapp.f  update.f
	$(CXX) $(CXXFLAGS) $<  newuoa/*.o -lgfortran -o $@

sync:
	cp -v condor-wrap.h gsl-nelder-mead-wrap.h Makefile  t.cc ../nf/lopti-ss/
