main: bigden.o biglag.o calfun.o main.o newuoa.o newuob.o trsapp.o update.o
	f77 -O -o newuoa bigden.o biglag.o calfun.o main.o newuoa.o newuob.o trsapp.o update.o
bigden.o: bigden.f
	f77 -O  -c bigden.f
biglag.o: biglag.f
	f77 -O  -c biglag.f
calfun.o: calfun.f
	f77 -O  -c calfun.f
main.o: main.f
	f77 -O  -c main.f
trsapp.o: trsapp.f
	f77 -O  -c trsapp.f
newuoa.o: newuoa.f
	f77 -O  -c newuoa.f
newuob.o: newuob.f
	f77 -O  -c newuob.f
update.o: update.f
	f77 -O  -c update.f

