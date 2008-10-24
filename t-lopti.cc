	
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

int main(int argc, char **argv) {
	
				typedef array<double,2,1>		array1_t;	
				typedef array<double,2,0>		array0_t;	

			array0_t		X0 = {{ -1.2, 1 }};
			//array0_t		X;
			of_rosenberg<array1_t>  of_rb1;
			of_rosenberg<array0_t>  of_rb0;

	{  ////  CONDOR  x PLAIN_FN  ROSENBERG
	condor_minimizer<array0_t>	mzr			(X0);	// X[0..N-1]
	mzr	.object_functOR		(new plain_fn<array0_t>(&plain_fn_rosenberg<array0_t>));
	mzr
		.rho_begin		(1)
		//.rho_end             (1e-3);
		.rho_end		(1e-10);
	mzr.argmin(); mzr.print();  }

	{  ////  CONDOR  x NAKED ROSENBERG
	condor_minimizer<array0_t>	mzr			(X0);	// X[0..N-1]
	mzr	.object_functOR		(&of_rb0)
		.rho_begin		(1)
		//.rho_end             (1e-3);
		.rho_end		(1e-10);
	mzr.argmin(); mzr.print();  }


	{  ////  CONDOR  logged rosenberg
	condor_minimizer<array0_t>	mzr			(X0);	// X[0..N-1]
	of_log<array0_t> ol  (&of_rb0, mzr);
	mzr	.object_functOR		(&ol)
	//	.object_functOR		(&of_log<array0_t>(&of_rb0,  mzr))
		.rho_begin		(1)
		//.rho_end             (1e-3);
		.rho_end		(1e-10);
	mzr.argmin(); mzr.print();  }	

	{  ////  NEWUOA :  2*N + 1  NAKED
	newuoa_minimizer<array1_t>	mzr			(*(array1_t*)&(X0));	// X[1..N]
	of_log<array1_t> ol  (&of_rb1, mzr);
	mzr	.object_functOR		(&of_rb1)
		.rho_begin		(1)
		.rho_end		(4e-11);
	mzr.argmin(); mzr.print(); }


	{  ////  NEWUOA :  2*N + 1 
	newuoa_minimizer<array1_t>	mzr			(*(array1_t*)&(X0));	// X[1..N]
	of_log<array1_t> ol  (&of_rb1, mzr);
	mzr	.object_functOR		(&ol)
		.rho_begin		(1)
		.rho_end		(4e-11);
	mzr.argmin(); mzr.print(); }


	{  ////  NEWUOA  (N+1)*(N+2)/2
	const int N=array1_t::sz;
	newuoa_minimizer<array1_t,(N+1)*(N+2)/2>	mzr			(*(array1_t*)&(X0));	// X[1..N]
	of_log<array1_t> ol  (&of_rb1, mzr);
	mzr	.object_functOR		(&ol)
		.rho_begin		(1)
		//.rho_end             (4e-4);
		.rho_end		(4e-11);
	mzr.argmin(); mzr.print(); }


/*
	{  ////  CONDOR (bad_scale_rosenbrock)
	const int FACTOR=100;
	of_bad_scale_rosenberg<array0_t, FACTOR>  of_bsrb0;
	array0_t X0 ={{-1.2,1*FACTOR}};
	condor_minimizer<array0_t>	mzr			(X0);	// X[0..N-1]
	mzr	.object_functION	(of_log<array0_t>(&of_bsrb0,  mzr))
		.object_functOR 	(&of_log<array0_t>(&of_bsrb0,  mzr))
		.rho_begin		(1)
		//.rho_end             (1e-3);
		.rho_end		(1e-10);
	mzr.argmin(); mzr.print(); }

								//array_t			R  = {{ 0.2, 0.2 }};
								//mzr.rescale		(R);
							
	{  ////  CONDOR (bad_scale_rosenbrock * rescale)
	const int FACTOR=100;
	of_bad_scale_rosenberg<array0_t, FACTOR>  of_bsrb0;
	array0_t X0 ={{-1.2,1*FACTOR}};
	condor_minimizer<array0_t>	mzr			(X0);	// X[0..N-1]
	mzr	.object_function	(of_rescale<array0_t>(&of_bsrb0,  mzr))
		.rho_begin		(1)
		//.rho_end             (1e-3);
		.rho_end		(1e-10);
	//mzr.argmin(); //mzr.rescale	(R); mzr.print(); }

 */
	{  ////  NELDER-MEAD
	array0_t		S  = {{ 0.6, 0.6 }};
	nelder_mead_minimizer<array0_t>	mzr			(X0);	// will ignore BEGIN index
	of_log<array0_t> ol  (&of_rb0, mzr);
	mzr	.object_functOR		(&ol)
		.step			(S)
		//.characteristic_size	(0.0002);
		.characteristic_size	(0.00000001);
	mzr.argmin(); mzr.print(); }

	return 0;
 }
