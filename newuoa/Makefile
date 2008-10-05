newuoa: bigden.o biglag.o calfun.o main.o newuoa.o newuob.o trsapp.o update.o
	gfortran -O -o newuoa bigden.o biglag.o calfun.o main.o newuoa.o newuob.o trsapp.o update.o
bigden.o: bigden.f
	gfortran -O  -c bigden.f
biglag.o: biglag.f
	gfortran -O  -c biglag.f
calfun.o: calfun.f
	gfortran -O  -c calfun.f
main.o: main.f
	gfortran -O  -c main.f
trsapp.o: trsapp.f
	gfortran -O  -c trsapp.f
newuoa.o: newuoa.f
	gfortran -O  -c newuoa.f
newuob.o: newuob.f
	gfortran -O  -c newuob.f
update.o: update.f
	gfortran -O  -c update.f

