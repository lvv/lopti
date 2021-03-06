	
	//////////////////////////////////////////////////////////////////// 
	#include <lvv/lvv.h>
	#include <lvv/math.h>
		using lvv::pow2;
	#include <lvv/array.h>
		using lvv::array;
		using lvv::matrix;

	#include	<lopti/object_function.h>
	#include	<lopti/condor-wrap.h>
	#include   	<lopti/newuoa-wrap.h>
	#include	<lopti/gsl-nelder-mead-wrap.h>
	#include	<lopti/hook-jeevs.h>

	using namespace lopti;

int main(int argc, char **argv) {
	
			///////////////////////////////////////////////// CONFIG
		
			#if  ! defined(NEWUOA) && ! defined(NM)  && !defined(CONDOR)  && !defined(PLAIN_FN)
				#define ALL
			#endif

			#ifdef  ALL
				#define CONDOR
				#define NEWUOA
				#define NM
				#define NAKED
				#define PLAIN_FN
			#endif


			#if  ! defined(FN) 
				#define FN rosenbrock
			#endif

			//#define STOP_AT_X_STEP 1e-3
			#define		STOP_AT_X_STEP 1e-15
				
			#if 	( FN == rosenbrock )
				#undef _N
				#define _N 2
			#endif

			#if  ! defined(_N) 
				#define _N 2
			#endif

			typedef 	array<double,_N,1>		V1;	
			typedef 	array<double,_N,0>		V0;	

			V0		_X0 = {{ -1.2, 1 }};
			#if 	(FN == chebyquad)
				for (int i=0; i <_N; i++)   _X0[i] = (i+1.)/(i+2.);
				const double RHO_BEGIN = 0.2* _X0[0];
				cout << "*** initial X: "<< _X0 << endl;
			#else
				#define 	RHO_BEGIN	0.1
			#endif

			
						/*
						cout << "@@@@@@@  " << rosenbrock<V1>().name() << endl; 
						cout << "@@@@@@@  " << chebyquad<V1>().name() << endl; 
						cout << "@@@@@@@  " << (trace<V1>(chebyquad<V1>())).name() << endl; 
						V1 RR;
						cout << "@@@@@@@  " << rescale<V1>
									(
										rosenbrock<V1>(),
										RR
									)
									.name() << endl; 
						cout << "@@@@@@@  " << rescale<V1>(rosenbrock<V1>(), (const V1)()) .name() << endl; 
						*/

	#ifdef CONDOR

		{  condor_minimizer<V0>	mzr;////  condor  LOGGED 
			mzr	.x0		(_X0);
			mzr	.objective	(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();	
		}

		#ifdef  NAKED
		{  condor_minimizer<V0>	mzr;	////  CONDOR  x NAKED ROSENBERG
			mzr	.objective	(  rosenbrock<V0>() );
			mzr	.x0		(_X0);
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}
		#endif

		#ifdef  RESCALE
		{  	V0 R = {{ 1, 0.0100 }};   V0 X = _X0; V0  X_opt;
		condor_minimizer<V0>	mzr;////  condor  logged RESCALED rosenbrock
			mzr	.objective	(xg_log<V0>  (rescale<V0>  ( trace<V0>(rosenbrock<V0>()), R),  mzr));
			mzr	.x0		(X/=R);	// X[0..N-1]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
		X_opt = mzr	.argmin();
		X_opt *= R;
			mzr	.print(); 

		cout << "optimum at: " << (X_opt *= R) << endl;
		}	
		#endif

	#endif


	#ifdef  NEWUOA

		/*
		{ newuoa_minimizer<V1>	mzr;		 ////  NEWUOA :  2*N + 1 
			mzr	.objective	(xg_log<V1>(FN<V1>(),mzr));
			mzr	.x0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin		();
			mzr	.print		();
		}*/

		{//  ANY V  (<float,2,0>)
		typedef		array<float,2,1>		V;
		V X0 = {{-1.2, 1}};
		X0 = _X0;

		newuoa_minimizer<V>	mzr;		 ////  NEWUOA :  2*N + 1 
			mzr	.objective	(xg_log<V>(FN<V>(),mzr));
			mzr	.x0		(X0);		// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin		();
			mzr	.print		();
		}



		#if _N < 15
		{	const int N=V1::sz;
		newuoa_minimizer<V0, (N+1)*(N+2)/2>	mzr; ////  NEWUOA  (N+1)*(N+2)/2
			mzr	.objective	(xg_log<V0>(FN<V0>(),mzr));
			mzr	.x0		(_X0);	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();
		}
		#endif

	#endif

	#ifdef  PLAIN_FN
	{  newuoa_minimizer<V0> mzr;
		mzr	.objective		( make_objective<V0> (&plain_fn_rosenbrock<V0>, "make_objective(plain_fn_rosenbrock)") );
		mzr	.x0			(_X0);	// X[0..N-1]	
		mzr	.rho_begin		(1);
		mzr	.rho_end		(STOP_AT_X_STEP);
		mzr	.argmin();
		mzr	.print(); 
		cout << "???????????????????????" << endl;
	}
	#endif


	#ifdef NM
		{	V0  S;   S = 0.1; //{{ 0.6, 0.6 }};
		gsl_nelder_mead_minimizer<V0>	mzr;	////  NELDER-MEAD
			mzr	.objective		(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.x0		(_X0);  // will ignore BEGIN index
			mzr	.step0		(S);
			//mzr	.characteristic_size	(0.0002);
			mzr	.characteristic_size	(STOP_AT_X_STEP);
			mzr.argmin();
			mzr.print();
		}
	#endif 

	#ifdef HJ
		{	V0  S;   S = 0.03; //{{ 0.6, 0.6 }};
		hook_jeevs_minimizer<V0>	mzr;	////  NELDER-MEAD
			mzr	.objective		(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.x0		(_X0);  // will ignore BEGIN index
			mzr	.step0		(S);
			//mzr	.characteristic_size	(0.0002);
			mzr	.tau		(1e-30);
			mzr.argmin();
			mzr.print();
		}
	#endif 

	cout << endl;
 }
