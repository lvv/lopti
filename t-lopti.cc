	
	#include	<lopti/lopti.h>
	//#include	<lopti/condor-wrap.h>
	//#include   	<lopti/newuoa-wrap.h>
	//#include	<lopti/gsl-nelder-mead-wrap.h>

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
	typename V::value_type&		// these refs make it work if X index start from 0 or 1
		x1 = *X.begin(), 
		x2 = *(X.begin()+1);
	//return 100*pow2(X[2]-pow2(X[1]))+pow2(1-X[1]); 
	assert(X.size() == 2);
	return 100*pow2(x2-pow2(x1))+pow2(1-x1); 
};

//#include <gzstream.h>
//#include <iostream>                                                                                                                                         
//#include <fstream>
	//using namespace std;
	//using std::ios;
	//using std::ofstream;
//#include <cstdlib>

			template<typename V>
struct	of_log  {
		typedef 	typename V::value_type		fp_t;
		typedef 	typename minimizer<V>::of_ptr_t		of_ptr_t;

		of_ptr_t	of_p;
		int		iter;
		//ogzstream	log_file;
		//ofstream	log_file;

	explicit	of_log		(of_ptr_t of, const char* log_name)		: of_p(of), iter(0) /*, log_file(log_name, ios::out)*/{};
	//explicit	~of_log		()						: {log_file.close()};

	fp_t	operator()	(V&  X)			{
		iter++;  
		fp_t   y = of_p(X); 
		//log_file << format("%d \t \%g \t %g \n") % iter %y %X; 
		return y;
	};

};

int main(int argc, char **argv) {
	
	typedef array<double,2,1>		array1_t;	
	typedef array<double,2,0>		array0_t;	
	array0_t		X0 = {{ -1.2, 1 }};
	array0_t		X;


	{
	newuoa_minimizer<array1_t>	mzr			(of_log<array1_t>(of_rosenberg<array1_t>, "ttt.gz"),  *(array1_t*)&(X0=X));	// X[1..N]
					mzr.rho_begin		(0.5);
					mzr.rho_end		(4e-4);
					mzr.argmin();
					mzr.print();
	}

	{
	condor_minimizer<array0_t>	mzr			(of_rosenberg<array0_t>, X0=X);	// X[0..N-1]
					mzr.rho_begin		(2);
					mzr.rho_end		(1e-3);
					mzr.argmin();
					mzr.print();
					//array_t			R  = {{ 0.2, 0.2 }};
					//mzr.rescale		(R);
	}

	{
	array0_t		S  = {{ 0.6, 0.6 }};
	nelder_mead_minimizer<array0_t>	mzr			(of_rosenberg<array0_t>,  X0=X);	// will ignore BEGIN index
					mzr.step		(S);
					mzr.characteristic_size	(0.0002);
					mzr.argmin();
					mzr.print();
	}

	return 0;
 }
