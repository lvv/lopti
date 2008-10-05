//
//     The Chebyquad test problem (Fletcher, 1965) for N = 2,4,6 and 8,
//     with NPT = 2N+1.
//

#include <lvv/lvv.h>
#include <lvv/array.h>
	using lvv::array;

extern "C"  void newuoa_ (int* N,  int* NPT,  double* X,  double* RHOBEG,  double* RHOEND,  int* IPRINT,  int* MAXFUN,  double* W);

int main() {

	array<double,10> X;
	array<double,10000> W;
	int IPRINT=2;
	int MAXFUN=5000;
	double RHOEND=1.0e-6;
	
	for (int N=2;  N <= 8;  N+=2)  {
                int NPT=2*N+1;

		for (int I=1; I<=N; I++) X[I]=I/double(N+1);

                double RHOBEG=0.2e0 * X[1];

                FMT("Results with N =%d and NPT =%d")  % N  % NPT;

              //CALL NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
              newuoa_ (&N, &NPT, (double*)&X, &RHOBEG, &RHOEND, &IPRINT, &MAXFUN, (double*) &W);
	}

return 0 ;}

