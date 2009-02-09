	
	/////////////////////////////////////////////////////////////////////////////// 
	#include	<lvv/lvv.h>
	#include	<lvv/math.h>
		using lvv::pow2;

	#include	<lvv/array.h>
		using lvv::array;
		using lvv::matrix;

	#include	<functional>
		using std::unary_function;
	#include	<boost/function.hpp>
		using boost::function;

	#include   	<lopti/newuoa-wrap.h>

	using namespace lopti;

int main(int argc, char **argv) {
	

			#if  ! defined(FN) 
				#define FN rosenberg
			#endif

			#define RHO_BEGIN 0.2
			#define RHO_END 1e-4
				
			#if  ! defined(_N) 
				#define _N 2
			#endif

			#if  ! defined(ITER) 
				#define ITER 1000
			#endif

			typedef 	array<double,_N,1>		V1;	
			typedef 	array<double,_N,0>		V0;	


		{ cout << "*****   ANY V  ****** \n";
		typedef		array<float,2,0>		V;
		V		X0 = {{ -1.2, -1 }};
		newuoa_minimizer<V>	mzr;		 ////  NEWUOA :  2*N + 1 
			//mzr	.loft		(xg_log<V>(FN<V>(),mzr));
			mzr	.loft		(trace<V>(FN<V>()));
			mzr	.x0		(X0);	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(RHO_END);
			mzr	.max_iter	(ITER);
			//mzr	.verbose	(true);
			V X_opt = mzr.argmin();
			mzr	.print		();
		}


	return 0;
 }
