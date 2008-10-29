	
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

	#include	<lopti/lopti.h>
		//#include	<lopti/condor-wrap.h>
		//#include   	<lopti/newuoa-wrap.h>
		//#include	<lopti/gsl-nelder-mead-wrap.h>
	#include <lopti/object_function.h>

	using namespace lopti;

int main(int argc, char **argv) {
	
				typedef 	array<double,_N,1>		V1;	
				typedef 	array<double,_N,0>		V0;	

			///////////////////////////////////////////////// CONFIG 
			#ifdef  ALL
				#define CONDOR
				#define NEWUOA
				#define NM
				#define NAKED
				#define PLAIN_FN
				#define  OLD_BAD_SCALE_ROSENBROCK
			#endif

			#if  ! defined(NEWUOA) && ! defined(NM)  && !defined(CONDOR)
				#define CONDOR
				#define NEWUOA
				#define NM
			#endif

			#if  ! defined(FN) 
				#define FN rosenberg
			#endif

			#if  ! defined(_N) 
				#define _N 2
			#endif
			
			//#define STOP_AT_X_STEP 1e-3
			#define STOP_AT_X_STEP 1e-10
				
			V0		_X0 = {{ -1.2, 1 }};

	#ifdef CONDOR

		{  condor_minimizer<V0>	mzr;////  condor  LOGGED chebyquad
			mzr	.X0		(_X0);
			mzr	.loft		(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();	
		}

		#ifdef  PLAIN_FN
		{  condor_minimizer<V0> mzr;////  CONDOR  x PLAIN_FN  ROSENBERG
			mzr	.loft			( plain_fn<V0> (&plain_fn_rosenberg<V0>, "plain_fn") );
			mzr	.X0			(_X0);	// X[0..N-1]	
			mzr	.rho_begin		(1);
			mzr	.rho_end		(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}
		#endif

		#ifdef  NAKED
		{  condor_minimizer<V0>	mzr;	////  CONDOR  x NAKED ROSENBERG
			mzr	.loft		(  rosenberg<V0>() );
			mzr	.X0		(_X0);
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}
		#endif


		#ifdef  RESCALE
		{  	V0 R = {{ 1, 0.0100 }};   V0 X = _X0; 
		condor_minimizer<V0>	mzr;////  condor  logged RESCALED rosenberg
			mzr	.loft		(xg_log<V0>  (rescale<V0>  (rosenberg<V0>(), R),  mzr));
			mzr	.X0		(X/=R);	// X[0..N-1]
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}	
		#endif

		
		#ifdef  OLD_BAD_SCALE_ROSENBROCK
		{ 	const int FACTOR=100;
			V0 _X0 ={{-1.2,1*FACTOR}};
		condor_minimizer<V0>	mzr;
			mzr	.X0		(_X0);	// X[0..N-1]
			mzr	.loft 		(xg_log<V0>  (bad_scale_rosenberg<V0, FACTOR>(),  mzr));
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();
		}
		#endif

	#endif


	#ifdef  NEWUOA

		{ newuoa_minimizer<V1>	mzr;		 ////  NEWUOA :  2*N + 1 
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.X0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin		();
			mzr	.print		();
		}


		{	const int N=V1::sz;
		newuoa_minimizer<V1, (N+1)*(N+2)/2>	mzr; ////  NEWUOA  (N+1)*(N+2)/2
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.X0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();
		}

	#endif


	#ifdef NM
		{	V0		S = {{0.6}}; //{{ 0.6, 0.6 }};
		nelder_mead_minimizer<V0>	mzr;	////  NELDER-MEAD
			mzr	.loft		(xg_log<V0>(FN<V0>(),  mzr));
			mzr	.X0		(_X0);  // will ignore BEGIN index
			mzr	.step		(S);
			//mzr	.characteristic_size	(0.0002);
			mzr	.characteristic_size	(STOP_AT_X_STEP);
			mzr.argmin();
			mzr.print();
		}
	#endif 

	return 0;
 }
