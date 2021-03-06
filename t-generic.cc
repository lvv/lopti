	
	/////////////////////////////////////////////////////////////////////////////// 
	#include	<lvv/lvv.h>
	#include	<lvv/math.h>
		using lvv::pow2;

	#include	<lvv/array.h>
		using lvv::array;
		using lvv::matrix;


	///////////////////////////////////////////////////////////////////  SELECT SOLVER
	#ifdef NM
		#include   	<lopti/gsl-nelder-mead-wrap.h>
	#endif

	#ifdef CONDOR
		#include   	<lopti/condor-wrap.h>
	#endif

	#ifdef HJ
		#include   	<lopti/hook-jeevs.h>
	#endif

	#ifdef NEWUOA
		#include   	<lopti/newuoa-wrap.h>
	#endif

	#ifndef MINIMIZER
		#include   	<lopti/newuoa-wrap.h>
		// #error "error:  no minimizer is selected"
	#endif
							

	///////////////////////////////////////////////////////////////////  SELECT OBJECTIVE
	
	#if  ! defined(FN) 
		#if  _N > 2 
			#define FN rosenbrock
		#else
			#define FN chebyquad
		#endif
	#endif

	#if  ! defined(_N) 
		#define _N 2
		//#warning 	"using default N=2"
	#endif



	using namespace lopti;

int main(int argc, char **argv) {
	
	typedef 	array<double,_N,0>		V;	
	V		X0 = {{ -1.2, 1 }};

	#if 	FN == chebyquad
		for (int i=0;  i < _N;  i++)   X0[i] = (i+1.)/(i+2.);
		//const double RHO_BEGIN = 0.2* X0[0];
		cout << "x0: "<< X0 << endl;
	#endif


	MINIMIZER<V>	mzr;
		mzr	.objective		(trace<V>(FN<V>()));
		mzr	.x0		(X0);	// X[1..N]
		#ifdef ITER
			mzr	.max_iter	(ITER);
		#else
			mzr	.max_iter	(500);
		#endif

		V X_opt = mzr.argmin();
		mzr	.print		();
		cout << endl;

 }
