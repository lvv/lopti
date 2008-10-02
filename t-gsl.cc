
	// for sleep:
	#include <unistd.h>
	    
	#include <stdio.h>
	#include <stdlib.h>
	#include <math.h>
	#include <string.h>

	#include <lvv/lvv.h>
	#include <lvv/math.h>
	using lvv::pow2;
	#include <lvv/array.h>
	using lvv::array;

	typedef array<double,2>		array_t;	

template<typename V>  typename V::value_type  of_rb(V X)   {  return 100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]);  };


int main(int argc, char **argv) {
	
	array_t		X0 = {{ -1.2, 1 }};
	array_t		R  = {{ 0.1, 0.1 }};

	/*{ 
		#include <lopti/condor_wrap.h>
		minimizer<array_t>	mzr(of_rb, X0);
		mzr.condor_rho_start	(1);
		mzr.condor_rho_end	(1e-7);
		mzr.rescale		(R);
		mzr.gnuplot_print	(true);
		array_t	Xmin = mzr.argmin();
		
		MSG("Result: Xmin%g   y=%g   iter=%d \n") %Xmin  %(mzr.ymin())  %(mzr.iter());
	}*/

	{ 
		#include <lopti/gsl-nelder-mead-wrap.h>
		minimizer<array_t>	mzr(of_rb, X0);
		mzr.condor_rho_start	(1);
		mzr.condor_rho_end	(1e-7);
		mzr.rescale		(R);
		mzr.gnuplot_print	(true);
		array_t	Xmin = mzr.argmin();
		
		MSG("Result: Xmin%g   y=%g   iter=%d \n") %Xmin  %(mzr.ymin())  %(mzr.iter());
	}

	return 0;
 }
