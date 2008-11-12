	
	/////////////////////////////////////////////////////////////////////////////// 
	#include <lvv/lvv.h>
	#include <lvv/math.h>
		using lvv::pow2;
	#include <lvv/array.h>
		using lvv::array;
		using lvv::matrix;

	// functional 
	#include <functional>
		//using std::binder1st;
		using std::unary_function;
	#include <boost/function.hpp>
		using boost::function;
	#include <boost/bind.hpp>
		using boost::bind;

	//#include	<lopti/object_function.h>
	//#include	<lopti/lopti.h>
	#include	<lopti/condor-wrap.h>
	#include   	<lopti/newuoa-wrap.h>
	#include	<lopti/gsl-nelder-mead-wrap.h>
	#include	<lopti/hook-jeevs.h>

	using namespace lopti;

int main(int argc, char **argv) {
	
			///////////////////////////////////////////////// CONFIG  <1>
			#ifdef  ALL
				#define CONDOR
				#define NEWUOA
				#define NM
				#define NAKED
				#define PLAIN_FN
			#endif

			#if  ! defined(NEWUOA) && ! defined(NM)  && !defined(CONDOR)
				#define CONDOR
				#define NEWUOA
				#define NM
				#define HJ
			#endif

			#if  ! defined(FN) 
				#define FN rosenberg
			#endif

			//#define STOP_AT_X_STEP 1e-3
			#define		STOP_AT_X_STEP 1e-15
				
			#if 	( FN == rosenberg )
		//		#undef _N
		//		#define _N 2
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
				cout << "x0: "<< _X0 << endl;
			#elif
				#define 	RHO_BEGIN	0.1
			#endif

			

	#ifdef CONDOR

		{  condor_minimizer<V0>	mzr;////  condor  LOGGED 
			mzr	.x0		(_X0);
			mzr	.loft		(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();	
		}

		#ifdef  PLAIN_FN
		{  condor_minimizer<V0> mzr;////  CONDOR  x PLAIN_FN  ROSENBERG		// TODO: why results deffrent from loft?
			//mzr	.loft			( make_loft<V0> (&plain_fn_rosenberg<V0>, "make_loft") );
			mzr	.loft			( make_loft<V0> (&bind (&plain_fn_rosenberg<V0>, _1) ));
			mzr	.x0			(_X0);	// X[0..N-1]	
			mzr	.rho_begin		(1);
			mzr	.rho_end		(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}
		#endif

		#ifdef  NAKED
		{  condor_minimizer<V0>	mzr;	////  CONDOR  x NAKED ROSENBERG
			mzr	.loft		(  rosenberg<V0>() );
			mzr	.x0		(_X0);
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}
		#endif

		#ifdef  RESCALE
		{  	V0 R = {{ 1, 0.0100 }};   V0 X = _X0; 
		condor_minimizer<V0>	mzr;////  condor  logged RESCALED rosenberg
			mzr	.loft		(xg_log<V0>  (rescale<V0>  ( trace<V0>(rosenberg<V0>()), R),  mzr));
			mzr	.x0		(X/=R);	// X[0..N-1]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}	
		#endif

	#endif


	#ifdef  NEWUOA

		{ newuoa_minimizer<V1>	mzr;		 ////  NEWUOA :  2*N + 1 
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.x0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin		();
			mzr	.print		();
		}


		#if _N < 15
		{	const int N=V1::sz;
		newuoa_minimizer<V1, (N+1)*(N+2)/2>	mzr; ////  NEWUOA  (N+1)*(N+2)/2
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.x0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();
		}
		#endif

	#endif


	#ifdef NM
		{	V0  S;   S.assign(0.1); //{{ 0.6, 0.6 }};
		gsl_nelder_mead_minimizer<V0>	mzr;	////  NELDER-MEAD
			mzr	.loft		(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.x0		(_X0);  // will ignore BEGIN index
			mzr	.step0		(S);
			//mzr	.characteristic_size	(0.0002);
			mzr	.characteristic_size	(STOP_AT_X_STEP);
			mzr.argmin();
			mzr.print();
		}
	#endif 

	#ifdef HJ
		{	V0  S;   S.assign(0.03); //{{ 0.6, 0.6 }};
		hook_jeevs_minimizer<V0>	mzr;	////  NELDER-MEAD
			mzr	.loft		(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.x0		(_X0);  // will ignore BEGIN index
			mzr	.step0		(S);
			//mzr	.characteristic_size	(0.0002);
			mzr	.tau		(1e-30);
			mzr.argmin();
			mzr.print();
		}
	#endif 

	return 0;
 }
