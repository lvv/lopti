	///////////////////////////////////////////////////////////////////////////////   LOPTI method selection
						// This selector only for this t-*.  
						// Application should select optimizer by including approprite *-wrap.h. 
						// Should be before array.h (so that gsl*.h  are before array.h)
	#ifndef LOPTI
	        //#define	LOPTI    NEWUOA
	        //#define	LOPTI    CONDOR
	        #define	LOPTI    NM
	#endif 
	
	#define         NEWUOA   	<lopti/newuoa-wrap.h>
	#define         NM              <lopti/gsl-nelder-mead-wrap.h>
	#define         CONDOR          <lopti/condor-wrap.h>
	        #include        LOPTI
	#undef          NM
	#undef          CONDOR
	#undef          NEWUOA

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

	//using namespace boost;
	//using namespace std;
	/////////////////////////////////////////////////////////////////////////////// OF

			template<typename V>  typename V::value_type
of_chebyquad		(V& X) 		{			// The Chebyquad test problem (Fletcher, 1965) 
	
	const int N = V::size();
	typedef  typename V::value_type fp_t;

	//array<array<fp_t,V::sz,1>, V::sz+1,1> Y;
	matrix <fp_t, V::sz, V::sz+1> Y;

     	for (int J=1; J<=N; J++)  {
		Y[1][J] = 1.0;
		Y[2][J] = 2.0*X[J]-1.0;
	}

     	for (int I=2; I<=N; I++) 
		for (int J=1; J<=N; J++)
			Y[I+1][J]=2.0*Y[2][J]*Y[I][J]-Y[I-1][J];

     	fp_t 	F  = 0.0;
     	int	NP = N+1;
     	int	IW = 1;

     	for (int I=1; I<=NP; I++)  {
		fp_t  SUM=0.0;
		for (int J=1; J<=N; J++) 	SUM += Y[I][J];
		SUM = SUM/N;
		if (IW > 0)  SUM += 1.0/(I*I-2*I);
		IW =-IW;
	   	F += SUM*SUM;
	}

	return F;
}

			template<typename V>  typename V::value_type
of_rosenberg		(V& X)   {  
	typename V::value_type&		// to make it work if X index start from 0 or 1
		x1 = *X.begin(), 
		x2 = *(X.begin()+1);
	//return 100*pow2(X[2]-pow2(X[1]))+pow2(1-X[1]); 
	assert(X.size() == 2);
	return 100*pow2(x2-pow2(x1))+pow2(1-x1); 
};



		//template<typename V>  typename V::value_type 
	//of_rb(V& X, void* var)   { return 100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]); };

	//double of_rb(array_t& X, void* var)   { return 100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]); };
	//double of_rb(array_t& X)   { return 100*pow2(X[1]-pow2(X[0]))+pow2(1-X[0]); };

int main(int argc, char **argv) {
	
	typedef array<double,2,1>		array_t;	
	typedef array<double,2,0>		array0_t;	
	array_t		X0 = {{ -1.2, 1 }};


			#if	defined(LOPTI_NEWUOA)
				newuoa_minimizer<array_t>	mzr			(of_rosenberg<array_t>,  X0);	// X[1..N]
								mzr.rho_begin		(0.5);
								mzr.rho_end		(1e-4);
								mzr.verbose		(true);
			#elif 	defined(LOPTI_CONDOR)
				condor_minimizer<array0_t>	mzr			(of_rosenberg<array0_t>, *(array0_t*)&X0);	// X[0..N-1]
								mzr.rho_begin		(2);
								mzr.rho_end		(1e-4);
								//array_t			R  = {{ 0.2, 0.2 }};
								//mzr.rescale		(R);
			#elif	defined(LOPTI_NM)
				array0_t		S  = {{ 0.6, 0.6 }};
				nelder_mead_minimizer<array0_t>	mzr			(of_rosenberg<array0_t>,  *(array0_t*)&X0);	// X[0..N-1]
								mzr.step		(S);
								mzr.gsl_characteristic_size(0.0001);
			#else
				assert(false && "macro LOPTI_<method>  not defined\n");
			#endif

				mzr.verbose		(true);
	array_t	Xmin =		*(array_t*) &(mzr.argmin());
	
	MSG("result: Xmin%g   y=%g   iter=%d  minimizer=%s\n") %Xmin  %(mzr.ymin())  %mzr.iter()  %mzr.name() ;

	return 0;
 }
