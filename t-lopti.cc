	
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

/////////////////////////////////////////////////////////////////////////////////////////  OF: CHEBYQUAD
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


/////////////////////////////////////////////////////////////////////////////////////////   OF_BASE

			template<typename V>
struct	of_base		{  
	of_base		(const char* s) :  		name_(s)   { X_opt_ = 0; };
	typedef  typename V::value_type	 fp_t;
	const char*	name_;
	V	X_opt_;
	int	iter_;

	virtual fp_t		operator()	(V& X) 		{ return 0; };
	virtual const char*	name		()	const	{ return  name_; };
	virtual int 		size		()	const	{ return  V::size(); };
	virtual int 		iter		()	const	{ return  iter_; };
	//virtual bool 		have_opt	()	const	{ return  false; };
		void 		opt		(const V& X_answ)	{ X_opt_ = X_answ; };
	virtual	fp_t 		opt_distance	(V& X)	const	{ return  distance_norm2(X_opt_, X); };
};


/////////////////////////////////////////////////////////////////////////////////////////  OF: ROSENBERG

			template<typename V>  
struct	of_rosenberg	: public of_base<V> { 
	of_rosenberg()	: of_base<V>("rosenberg") 	{ V const  X_answ = {{ 1.0, 1.0 }};   of_base<V>::opt(X_answ); };

			virtual  typename V::value_type
	operator() 	(V& X)   {  
		typename V::value_type&		// these refs make it work if X index start from 0 or 1
			x1 = *X.begin(), 
			x2 = *(X.begin()+1);
		return 100*pow2(x2-pow2(x1))+pow2(1-x1); 
	};
};

/////////////////////////////////////////////////////////////////////////////////////////  OF WRAP:  OF_LOG
//#include <gzstream.h>
#include <iostream>                                                                                                                                         
#include <fstream>
	using std::ios;
	using std::ofstream;

			template<typename V>
struct	of_log  {
		typedef 	typename V::value_type			fp_t;
		typedef 	typename minimizer<V>::of_ptr_t		of_ptr_t;

		of_base<V>*	of_base_p;
		of_ptr_t	of_p;
		int		iter;
		ofstream&	log_file;

	explicit	of_log		(of_base<V>* of, ofstream& lf)		: of_p(*of), of_base_p((of_base<V>*)of),  iter(0), log_file(lf)  { assert(log_file.good());};
	//explicit	of_log		(of_ptr_t of, ofstream& lf)		: of_p(of), of_base_p(&of),  iter(0), log_file(lf)  { assert(log_file.good());};

	fp_t	operator()	(V&  X)			{
		iter++;  
		//fp_t   y = of_p(X); 
		fp_t   y = of_base_p->operator()(X); 
								assert(log_file.good());
		//log_file << format("%d \t %g \t %g \n") % iter %y %X; 
		log_file << format("%d \t %g \t %g %g  %g\n") % iter %y  %of_base_p->opt_distance(X) %X %of_base_p->X_opt_;
		return y;
	};
};

int main(int argc, char **argv) {
	
	typedef array<double,2,1>		array1_t;	
	typedef array<double,2,0>		array0_t;	
	array0_t		X0 = {{ -1.2, 1 }};
	array0_t		X;
	ofstream		log_file("log");
	of_rosenberg<array1_t>  of_rb;

	{
	newuoa_minimizer<array1_t>	mzr			(of_log<array1_t>(&of_rb, log_file),  *(array1_t*)&(X=X0));	// X[1..N]
	//newuoa_minimizer<array1_t>	mzr			(of_log<array1_t>(of_rosenberg<array1_t>(), log_file),  *(array1_t*)&(X=X0));	// X[1..N]
					mzr.rho_begin		(0.5);
					mzr.rho_end		(4e-4);
					mzr.argmin();
					mzr.print();
	}

	{
	condor_minimizer<array0_t>	mzr			(of_rosenberg<array0_t>(), X=X0);	// X[0..N-1]
					mzr.rho_begin		(2);
					mzr.rho_end		(1e-3);
					mzr.argmin();
					mzr.print();
					//array_t			R  = {{ 0.2, 0.2 }};
					//mzr.rescale		(R);
	}

	{
	array0_t		S  = {{ 0.6, 0.6 }};
	nelder_mead_minimizer<array0_t>	mzr			(of_rosenberg<array0_t>(),  X=X0);	// will ignore BEGIN index
					mzr.step		(S);
					mzr.characteristic_size	(0.0002);
					mzr.argmin();
					mzr.print();
	}

	return 0;
 }
