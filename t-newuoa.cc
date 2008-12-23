	
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

			//#define STOP_AT_X_STEP 1e-3
			#define		STOP_AT_X_STEP 1e-15
				
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
			#else
				#define 	RHO_BEGIN	0.1
			#endif



		{ newuoa_minimizer<V1>	mzr;		 ////  NEWUOA :  2*N + 1 
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.x0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			//cout << mzr	.argmin		() << endl;
			mzr	.print		();
		}


		#if _N < 15
		{	const int N=V1::sz;
		newuoa_minimizer<V1, (N+1)*(N+2)/2>	mzr; ////  NEWUOA  (N+1)*(N+2)/2
			mzr	.loft		(xg_log<V1>(FN<V1>(),mzr));
			mzr	.x0		(*(V1*)&(_X0));	// X[1..N]
			mzr	.rho_begin	(RHO_BEGIN);
			mzr	.rho_end	(STOP_AT_X_STEP);
			//cout << mzr	.argmin() << endl;
			mzr	.print();
		}
		#endif


	return 0;
 }
