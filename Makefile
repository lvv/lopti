
include  ../lvv/include.mk

.DEFAULT_GOAL := t-condor
	
LDFLAGS += -lgsl -lgslcblas -L /usr/local/lib -lcondor
CXXFLAGS += -I /home/lvv/NF/

t-condor: CXXFLAGS +=   -I ..
t-condor: LDFLAGS  +=   -L /usr/local/lib/ -lcondor  -lm 

t-condor: t.cc condor-wrap.h
	$(CXX) $(CXXFLAGS) -DOPTI=CONDOR  $< -o $@ $(LDFLAGS)

t-nm: t.cc  gsl-nelder-mead-wrap.h
	$(CXX) $(CXXFLAGS) -DOPTI=NM  $< -o $@ $(LDFLAGS)

install:
	cd /usr/local &&  git checkout lvvlib 
	mkdir -p /usr/local/include/lvvlib/
	cp *.h   /usr/local/include/lvvlib/
	cd /usr/local &&  git add /usr/local/include/lvvlib
	cd /usr/local &&  git commit -a -m up
	cd /usr/local &&  git checkout installed
	cd /usr/local &&  git merge lvvlib


sync:
	cp -v condor-wrap.h gsl-nelder-mead-wrap.h Makefile  t.cc ../nf/lopti-ss/
