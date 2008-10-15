
	#include "condor-wrap.h"

	#include <lvv/math.h>
		using lvv::pow2;


			template<typename V>  typename V::value_type 
	of_rb(V& X, void* var)   {  return 100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]);  };


int main(int argc, char **argv) {
	
	typedef lvv::array<double,2>		array_t;	

	array_t		X0 = {{ -1.2, 1 }};
	array_t		R  = {{ 2, 2 }};

	minimizer<array_t>	mzr(of_rb, X0);
	mzr.condor_rho_start	(1);
	mzr.condor_rho_end	(1e-4);
	mzr.rescale		(R);
	mzr.verbose		(true);
	array_t	Xmin = mzr.argmin();
	
	MSG("# Result: Xmin%.10g   y=%.10g   iter=%d \n") %Xmin  %(mzr.ymin())  %(mzr.iter());

	return 0;
 }
