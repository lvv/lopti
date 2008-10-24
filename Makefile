
include  ../lvv/include.mk
	
LDFLAGS += -lgsl -lgslcblas -L /usr/local/lib -lcondor newuoa/*.o -lgfortran 
CXXFLAGS += -I /home/lvv/NF/
XGRAPHIC = xgraphic  -mark -markcol=-1  -g2 -logy -leg -legpos=3  -legsiz=1 -legtyp=2 -titgen="Convergance speed for dirivative-free algorithms" -titx="Objective function evaluation count" -tity="Distance to optimum: log10(|X-X_opt|)" 

.DEFAULT_GOAL := t-lopti-r

t-lopti-xg: t-lopti log/condor
	$<
	cd log; $(XGRAPHIC) *

log/condor: t-lopti 
	rm -f log/*

t-lopti-r: log/condor

t-lopti.ps: t-lopti log/condor
	cd log; $(XGRAPHIC) -pscoltex=../$@ *
	gsview $@

t-lopti: t-lopti.cc *.cc *.h

t-condor: CXXFLAGS +=   -I ..
t-condor: LDFLAGS  +=   -L /usr/local/lib/ -lcondor  -lm 

t-condor: t-condor.cc condor-wrap.h
	$(CXX) $(CXXFLAGS) -DOPTI=CONDOR  $< -o $@ $(LDFLAGS)

t-nm: t.cc  gsl-nelder-mead-wrap.h
	$(CXX) $(CXXFLAGS) -DOPTI=NM  $< -o $@ $(LDFLAGS)

git-install:
	cd /usr/local &&  git checkout lvvlib 
	mkdir -p /usr/local/include/lvvlib/
	cp *.h   /usr/local/include/lvvlib/
	cd /usr/local &&  git add /usr/local/include/lvvlib
	cd /usr/local &&  git commit -a -m up
	cd /usr/local &&  git checkout installed
	cd /usr/local &&  git merge lvvlib


#FCFLAGS = -frange-check -fbounds-check -O0 -ggdb3

t-newuoa: newuoa/t-newuoa.cc newuoa-wrap.h  lopti.h
	#	rm -f *.o
	#$(FC)  $(FCFLAGS) -c bigden.f  biglag.f  calfun.f   trsapp.f  update.f
	$(CXX) $(CXXFLAGS) $<  newuoa/*.o -lgfortran -o $@

sync:
	cp -v condor-wrap.h gsl-nelder-mead-wrap.h Makefile  t.cc ../nf/lopti-ss/
