
include  /home/lvv/NF/lvv/include.mk

CXXFLAGS += -ggdb3 -I /home/lvv/NF/
FC = gfortran
FCFLAGS = -frange-check -fbounds-check -O0 -ggdb3

t-newuoa: t-newuoa.cc newuoa-wrap.h  ../lopti.h bigden.f  biglag.f  calfun.f  main.f  newuoa.f  newuob.f  trsapp.f  update.f
	rm -f *.o
	$(FC)  $(FCFLAGS) -c bigden.f  biglag.f  calfun.f   trsapp.f  update.f
	$(CXX) $(CXXFLAGS) $< *.o -lgfortran -o $@

olduoa: olduoa.cc bigden.f  biglag.f  calfun.f  main.f  newuoa.f  newuob.f  trsapp.f  update.f
	rm -f *.o
	$(FC)  $(FCFLAGS) -c bigden.f  biglag.f  calfun.f   trsapp.f  update.f
	$(CXX) $(CXXFLAGS) $< *.o -lgfortran -o $@

newuoa: newuoa.cc bigden.f  biglag.f  calfun.f  main.f  newuoa.f  newuob.f  trsapp.f  update.f
	rm -f *.o
	$(FC)  $(FCFLAGS) -c bigden.f  biglag.f  calfun.f   trsapp.f  update.f
	$(CXX) $(CXXFLAGS) $< *.o -lgfortran -o $@

#	./$@

f-main: main.f  bigden.f  biglag.f  calfun.f  newuoa.f  newuob.f  trsapp.f  update.f
	rm -f *.o 
	$(FC)  $(FCFLAGS) bigden.f  biglag.f  calfun.f  main.f  newuoa.f  newuob.f  trsapp.f  update.f    -lgfortran  -o $@

c-main: main.c *.f
	rm -f *.o 
	gfortran -c *.f
	rm -f main.o 
	f2c main.f 
	$(CC) $< *.o -lf2c -lgfortran  -o $@

cc-main: main.cc *.f
	rm -f *.o 
	gfortran -c *.f
	rm -f main.o 
	$(CXX) $< *.o -lf2c -lgfortran  -o $@

main.c: main.f
	f2c -a -f main.f

#main.cc: main.f
#	f2c -a -f main.f -C++
#	mv main.c main.cc
