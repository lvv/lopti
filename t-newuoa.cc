	
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
			/*
			#if 	(FN == chebyquad)
				for (int i=0; i <_N; i++)   _X0[i] = (i+1.)/(i+2.);
				const double RHO_BEGIN = 0.2* _X0[0];
				cout << "x0: "<< _X0 << endl;
			#else
				#define 	RHO_BEGIN	0.1
			#endif
			*/



		{ newuoa_minimizer<V1>	mzr;		 ////  NEWUOA :  2*N + 1 
			V1		_X0 = {{ -1.2, 1 }};
			//mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.loft		(FN<V1>());
			mzr	.x0		(_X0);	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(RHO_END);
			mzr	.max_iter	(ITER);
			mzr	.verbose     (true);
			V1 X_opt = mzr.argmin();
			//cout << mzr	.argmin		() << endl;
			mzr	.print		();
		}

		{ newuoa_minimizer<V1>	mzr;		 ////  NEWUOA :  2*N + 1 
			V0		_X0 = {{ -1.2, 1 }};
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.x0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(RHO_END);
			mzr	.verbose     (true);
			mzr	.max_iter	(ITER);
			cout << mzr	.argmin		() << endl;
			mzr	.print		();
		}


		#if _N < 15
		{	const int N=V1::sz;
		newuoa_minimizer<V1, (N+1)*(N+2)/2>	mzr; ////  NEWUOA  (N+1)*(N+2)/2
			V0		_X0 = {{ -1.2, 1 }};
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.x0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(RHO_END);
			mzr	.verbose     (true);
			mzr	.max_iter	(ITER);
			cout << mzr	.argmin() << endl;
			mzr	.print();
		}
		#endif


	return 0;
 }
