	
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
	//#include "boost/noncopyable.hpp"
	//	using boost::noncopyable;
	#include <boost/bind.hpp>
		using boost::bind;

	#include	<lopti/lopti.h>
		//#include	<lopti/condor-wrap.h>
		//#include   	<lopti/newuoa-wrap.h>
		//#include	<lopti/gsl-nelder-mead-wrap.h>
	#include <lopti/object_function.h>

	using namespace lopti;

int main(int argc, char **argv) {
	
				typedef 	array<double,2,1>		array1_t;	
				typedef 	array<double,2,0>		array0_t;	

			//array0_t		X;
			of_rosenberg<array1_t>  of_rb1;
			of_rosenberg<array0_t>  of_rb0;

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
			#endif
			
			//#define STOP_AT_X_STEP 1e-3
			#define STOP_AT_X_STEP 1e-10
				
			array0_t		_X0 = {{ -1.2, 1 }};

	#ifdef CONDOR

		#ifdef  PLAIN_FN
		{  condor_minimizer<array0_t> mzr;////  CONDOR  x PLAIN_FN  ROSENBERG
			mzr	.loft			( plain_fn<array0_t> (&plain_fn_rosenberg<array0_t>, "plain_fn") );
			mzr	.X0			(_X0);	// X[0..N-1]	
			mzr	.rho_begin		(1);
			mzr	.rho_end		(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}
		#endif

		#ifdef  NAKED
		{  condor_minimizer<array0_t>	mzr;	////  CONDOR  x NAKED ROSENBERG
			mzr	.loft		(  of_rosenberg<array0_t>() );
			mzr	.X0		(_X0);
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}
		#endif


		{  condor_minimizer<array0_t>	mzr;////  condor  LOGGED rosenberg
			mzr	.X0		(_X0);
			mzr	.loft		(of_log<array0_t>(of_rosenberg<array0_t>(),  mzr));
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();	
		}	


		{  	array0_t R = {{1,1}};
		condor_minimizer<array0_t>	mzr;////  condor  logged RESCALED rosenberg
			mzr	.loft		(of_log<array0_t>  (rescale<array0_t>  (of_rosenberg<array0_t>(), R),  mzr));
			mzr	.X0		(_X0);	// X[0..N-1]
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print(); 
		}	

		
		#ifdef  OLD_BAD_SCALE_ROSENBROCK
		{ 	const int FACTOR=100;
			array0_t _X0 ={{-1.2,1*FACTOR}};
		condor_minimizer<array0_t>	mzr;
			mzr	.X0		(_X0);	// X[0..N-1]
			mzr	.loft 		(of_log<array0_t>  (of_bad_scale_rosenberg<array0_t, FACTOR>(),  mzr));
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();
		}
		#endif

	#endif


	#ifdef  NEWUOA

		{ newuoa_minimizer<array1_t>	mzr;		 ////  NEWUOA :  2*N + 1 
			mzr	.loft		(of_log<array1_t>(of_rosenberg<array1_t>(),mzr));
			mzr	.X0		(*(array1_t*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin		();
			mzr	.print		();
		}


		{	const int N=array1_t::sz;
		newuoa_minimizer<array1_t, (N+1)*(N+2)/2>	mzr; ////  NEWUOA  (N+1)*(N+2)/2
			mzr	.loft		(of_log<array1_t>(of_rosenberg<array1_t>(),mzr));
			mzr	.X0		(*(array1_t*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(1);
			mzr	.rho_end	(STOP_AT_X_STEP);
			mzr	.argmin();
			mzr	.print();
		}

	#endif


	#ifdef NM
		{	array0_t		S  = {{ 0.6, 0.6 }};
		nelder_mead_minimizer<array0_t>	mzr;	////  NELDER-MEAD
			mzr	.loft		(of_log<array0_t>(of_rosenberg<array0_t>(),  mzr));
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
