CXXFLAGS = -I /home/lvv/NF/

main: main.cc *.f
	gfortran -c *.f
	rm -f main.o t.o
	$(CXX) $(CXXFLAGS) $< *.o -lgfortran -o $@

#	$(CXX) $(CXXFLAGS) $< *.o -lf2c -lgfortran  -o $@

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
