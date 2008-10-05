



	///////////////////////////////////////////////////////////////////////////////   OPTI selection
						// should be before array.h
	#ifndef OPTI
	        #define OPTI    NM
	#endif 
	
	#define         NM              <lopti/gsl-nelder-mead-wrap.h>
	#define         CONDOR          <lopti/condor-wrap.h>
	        #include        OPTI
	#undef          NM
	#undef          CONDOR

	/////////////////////////////////////////////////////////////////////////////// 
	#include <lvv/lvv.h>
	#include <lvv/math.h>
		using lvv::pow2;
	#include <lvv/array.h>
		using lvv::array;

	/////////////////////////////////////////////////////////////////////////////// OF

	typedef array<double,2>		array_t;	

		template<typename V>  typename V::value_type 
	of_rb(V& X, void* var)   { return 100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]); };

int main(int argc, char **argv) {
	
	array_t		X0 = {{ -1.2, 1 }};

	minimizer<array_t>	mzr			(of_rb, X0);
			#ifdef OPTI_CONDOR
				mzr.condor_rho_start	(2);
				mzr.condor_rho_end	(1e-4);
				array_t			R  = {{ 0.2, 0.2 }};
				mzr.rescale		(R);
			#else
				array_t		S  = {{ 0.5, 0.5 }};
				mzr.step		(S);
				mzr.gsl_characteristic_size(0.00001);
			#endif
				mzr.verbose		(true);
	array_t	Xmin =		mzr.argmin		();
	
	MSG("Result: Xmin%g   y=%g   iter=%d \n") %Xmin  %(mzr.ymin())  %(mzr.iter());

	return 0;
 }
