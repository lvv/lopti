	
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
	//#include "boost/noncopyable.hpp"
	//	using boost::noncopyable;
	#include <boost/bind.hpp>
		using boost::bind;
	#include <boost/shared_ptr.hpp>
		using boost::shared_ptr;


	//using namespace boost;
	//using namespace std;

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

	virtual fp_t		operator()	(V& X) 		{ assert(false); return X[0]; };
	virtual const char*	name		()	const	{ return  name_; };
	virtual int 		size		()	const	{ return  V::size(); };
	virtual int 		iter		()	const	{ return  iter_; };
		void 		opt		(const V& X_answ)	{ X_opt_ = X_answ; };
	virtual	fp_t 		opt_distance	(V& X)	const	{ return  distance_norm2(X_opt_, X); };
	virtual void		reset		()		{ iter_ = 0;};
};


/////////////////////////////////////////////////////////////////////////////////////////  OF: ROSENBERG

			template<typename V>  
struct	of_rosenberg	: public of_base<V> { 
	of_rosenberg()	: of_base<V>("rosenberg") 	{ V const  X_answ = {{ 1.0, 1.0 }};   of_base<V>::opt(X_answ); };
	static	const int B = 		V::ibg;

			virtual  typename V::value_type
	operator() 	(V& X)   {  
		of_base<V>::iter_++;
		return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); 
	};
};

/////////////////////////////////////////////////////////////////////////////////////////  OF WRAP:  OF_LOG

//#include <gzstream.h>
#include <iostream>                                                                                                                                         
#include <fstream>
	using std::ios;
	using std::ofstream;

			template<typename V>
//class	of_log  : private boost::noncopyable { public:
class	of_log  { public:
		typedef 		typename V::value_type			fp_t;
		typedef 		typename minimizer<V>::of_ptr_t		of_ptr_t;

		of_base<V>*		of_base_p;
		of_ptr_t		of_p;
		int			iter;
		shared_ptr<ofstream> log_file;

	//explicit	of_log		(of_base<V>* of, ofstream& lf, const char* title)	
	explicit	of_log		(of_base<V>* of, const char* title)	
	:	of_p(*of),
		of_base_p((of_base<V>*)of), 
		iter(0),
		log_file(new ofstream(("log/"+string(title)).c_str()))
	{
									assert (log_file->good());
		of_base_p->reset();
		// SHOULD BE set grid; plot "1" using 1:(log10($2)) with lines,"2" using 1:(log10($2)) w l,"3" using 1:(log10($2)) w l

		// plot:  opt_dist (iter)
		*log_file << format("# :gnuplot+: \"pipe\" using 1:(log10($2)) with lines title \"%s(%s)\",  \n")  %title % of_base_p->name();

		// splot:  path 
		#undef	 GP_F
		#define  GP_F "set view 0,0,1.7;  splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		*log_file << format(
			"# :gnuplot2: set font \"arial,6\"; set dgrid3d;  set key off;"
			"set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 4:5:2:1 with labels title \"%s(%s)\" ; \n" 
		)  %title % of_base_p->name();
		 
		// col header
		*log_file << "# iter \t y=f(X) \t |y-opt| \t X\n";
	};

	//~of_log		() { log_file.close(); };


	fp_t	operator()	(V&  X)			{
		iter++;  
		//fp_t   y = of_p(X); 
		fp_t   y = (*of_base_p)(X); 
								assert(log_file->good());
		// 1 - iter no(gp: splot label);  2 - hight (y);  3 - |X-opt| ignored;  4,5 - X coord;  
		*log_file << format("%d \t %g \t %g \t %g \n") % iter  %of_base_p->opt_distance(X) %y %X; 
		return y;
	};
};

int main(int argc, char **argv) {
	
	typedef array<double,2,1>		array1_t;	
	typedef array<double,2,0>		array0_t;	
	array0_t		X0 = {{ -1.2, 1 }};
	//array0_t		X;
	of_rosenberg<array1_t>  of_rb1;
	of_rosenberg<array0_t>  of_rb0;


	////  NEWUOA
	newuoa_minimizer<array1_t>	no_mzr			(*(array1_t*)&(X0));	// X[1..N]
	no_mzr	.object_function	(of_log<array1_t>(&of_rb1, "newuoa"))
		.rho_begin		(2)
		.rho_end		(4e-11);
	no_mzr.argmin();
	no_mzr.print();		cout << " of_iter=" << of_rb1.iter() << endl;

	////  CONDOR
	condor_minimizer<array0_t>	c_mzr			(X0);	// X[0..N-1]
	c_mzr	.object_function	(of_log<array0_t>(&of_rb0,  "condor"))
		.rho_begin		(2)
		.rho_end		(1e-10);
	c_mzr.argmin();
	c_mzr.print();		cout << " of_iter=" << of_rb0.iter() << endl;
	//array_t			R  = {{ 0.2, 0.2 }};
	//mzr.rescale		(R);

	////  NELDER-MEAD
	array0_t		S  = {{ 0.6, 0.6 }};
	nelder_mead_minimizer<array0_t>	nm_mzr			(X0);	// will ignore BEGIN index
	nm_mzr	.object_function	(of_log<array0_t>(&of_rb0, "nelder_mead"))
		.step			(S)
		//.characteristic_size	(0.0002);
		.characteristic_size	(0.00000001);
	nm_mzr.argmin();
	nm_mzr.print();		cout << " of_iter=" << of_rb0.iter() << endl;
	return 0;
 }
