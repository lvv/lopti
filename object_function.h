
	#include <boost/shared_ptr.hpp>
		using boost::shared_ptr; // used as ofstream ptr

	//#include <gzstream.h>
	#include <iostream>                                                                                                                                         
	#include <fstream>
		using std::ios;
		using std::ofstream;


namespace lopti {
/////////////////////////////////////////////////////////////////////////////////////////  OF: CHEBYQUAD
			template<typename V>  typename V::value_type
of_chebyquad		(V& X) 		{			// The Chebyquad test problem (Fletcher, 1965) 
	
	const int N = V::size();
	typedef  typename V::value_type fp_t;

	matrix <fp_t, V::sz, V::sz+1>		Y;

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


/////////////////////////////////////////////////////////////////////////////////////////  ADAPTER: PLAIN_FN

			template<typename V>
struct	plain_fn : public loft_base<V>		{  	
					typedef		typename V::value_type		fp_t;
				function<fp_t(V&)>	of;
	explicit		plain_fn	(function<fp_t(V&)> _of, const string& _name)	 : loft_base<V>(_name) , of(_of)   {};
	virtual plain_fn<V>&	clone		()	const		{  cout << "plain_fn.clone\n";  return  *new plain_fn<V>(*this); }
	virtual fp_t		operator()	(V&  X)	  { assert(!of.empty() && ">> NOT DEFINED OBJ FUNC <<");  this->iter_++;   fp_t y = (of)(X); return y; }
 };


/////////////////////////////////////////////////////////////////////////////////////////  OF: ROSENBERG
			template<typename V>  
struct	of_rosenberg	: loft_base<V> { 
					static const int B = 		V::ibg;
					typedef		typename V::value_type		fp_t;
					typedef		of_rosenberg<V>	this_t;
	explicit 			of_rosenberg	()	: loft_base<V>("rosenberg") 	{ V const  X_answ = {{ 1.0, 1.0 }};   loft_base<V>::opt(X_answ); };
	virtual	this_t&			clone		()	const		{  cout << "rb.clone\n"; return  *new this_t(*this); }
	virtual  fp_t			operator() 	(V& X)  {  loft_base<V>::iter_++;      return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };
 };


template<typename V>  typename V::value_type plain_fn_rosenberg(V& X) { static const int B = V::ibg;   return  100 * pow2(X[1+B]-pow2(X[0+B])) + pow2(1-X[0+B]); };


/////////////////////////////////////////////////////////////////////////////////////////  OF: BAD SCALE ROSENBERG
			template<typename V, int FACTOR>  
struct	of_bad_scale_rosenberg	 : loft_base<V> { 
					typedef		typename V::value_type		fp_t;
					static	const int B = V::ibg;
					typedef			      loft_base<V>*		loft_p_t;
					typedef 		const loft_base<V>&		loft_cref_t;
					typedef			of_bad_scale_rosenberg<V,FACTOR> this_t;
	explicit 		of_bad_scale_rosenberg	() : loft_base<V>("bad_scale_rosenberg") 	{ V const  X_answ = {{ 1.0, 1.0*FACTOR }};   loft_base<V>::opt(X_answ); };
	virtual	this_t&		clone			()	const		{  cout << "rescale::clone()\n"; return  *new this_t(*this); }

			virtual  fp_t
	operator() 	(V& X)   {  
		loft_base<V>::iter_++;
		return  100 * pow2(X[1+B]/FACTOR-pow2(X[0+B])) + pow2(1-X[0+B]); 
	};
 };

/////////////////////////////////////////////////////////////////////////////////////////  OF: RESCALER


			template<typename V>
class	rescale :  public loft_base<V>  { public:
					typedef 		typename V::value_type			fp_t;
					typedef			      loft_base<V>*		loft_p_t;
					typedef 		const loft_base<V>&		loft_cref_t;
					typedef			rescale<V>			this_t;
			V R;
					explicit
	rescale		(const loft_base<V>& loft_v, V& _R)  : 	R(_R)   { loft_base<V>::loft(loft_v); };
	virtual	this_t&			clone		()	const		{  cout << "rescale::clone()\n"; return  *new this_t(*this); }

	virtual const string	name		()	const	{ return  *new string(this->name_ + " rescaled"); };  // FIXME leak
	//virtual const string	name		()		{ return  this->name_ = *new string(this->wrapped_loft_v->name_ +  " rescaled"); };
	//virtual const string	name		()		{ return  this->wrapped_loft_v->name_ +  " rescaled"; };

	fp_t	operator()	(V&  X)			{
		loft_base<V>::iter_++;  
		V XR = X; // so that here are no side effects of modified X
		XR *= R;
		fp_t   y = (*this->wrapped_loft_v)(XR); 
		return y;
	};
	
	virtual	fp_t 		opt_distance	(V& X)	const	{ V XR = X; XR*=R; return  distance_norm2(this->X_opt_, XR); };
 };
/////////////////////////////////////////////////////////////////////////////////////////  OF WRAP:  OF_LOG

			template<typename V>
class	of_log : public loft_base<V> { public:
					// takes ob_base-type functor, and calling minimizer (to extract its name)
					// and  logs all evals to file, which name constructed from minimizer-name and functin-name
					typedef 		typename V::value_type		fp_t;
					typedef			      loft_base<V>*		loft_p_t;
					typedef 		const loft_base<V>&		loft_cref_t;
					typedef			of_log<V>			this_t;

				shared_ptr<ofstream>	log_file;

	virtual	this_t&			clone		()	const		{  cout << "of_log::clone()\n"; return  *new this_t(*this); }


					explicit
	of_log	 (loft_cref_t _loft_v, minimizer<V>& mzr)	
	//:	loft_base<V>(loft_v),  		// init parent 
	//:	log_file(new ofstream(("log/" +  (&mzr)->name() + "(" + loft_base<V>::name() + ")" ).c_str()))
	{
		loft_base<V>::loft(_loft_v),  		// init parent 
		cout << "filename:  log/" +  (&mzr)->name() + "(" + loft_base<V>::name() + ")" << endl;
		log_file = shared_ptr<ofstream>(new ofstream(("log/" +  (&mzr)->name() + "(" + loft_base<V>::name() + ")" ).c_str()));
									assert (log_file->good());
		//this->wrapped_loft_v->reset();
		// SHOULD BE set grid; plot "1" using 1:(log10($2)) with lines,"2" using 1:(log10($2)) w l,"3" using 1:(log10($2)) w l

		// plot:  opt_dist (iter)
		//*log_file << format("# :gnuplot+: \"pipe\" using 1:(log10($2)) with lines title \"%s(%s)\",  \n")  %title % wrapped_loft_v->name();

		// splot:  path 
		/*
		#undef	 GP_F
		#define  GP_F "set view 0,0,1.7;  splot [-2:1.5][-0.5:2] log(100 * (y - x*x)**2 + (1 - x)**2),  "
		*log_file << format(
			"# :gnuplot2: set font \"arial,6\"; set dgrid3d;  set key off;"
			"set contour surface;  set cntrparam levels 20;  set isosample 40;"
			GP_F "\"pipe\" using 4:5:2:1 with labels title \"%s(%s)\" ; \n" 
		)  %title % wrapped_loft_v->name();
		 
		// col header
		*log_file << "# iter \t y=f(X) \t |y-opt| \t X\n";
		*/
	};

	~of_log		() { if(log_file.unique())   log_file->close(); };  // can not close file at DTOR.  On copy ctor,  of_log obj will be destoyed several times
	//void close		() { log_file->close(); };  // can not close file at DTOR.  On copy ctor,  of_log obj will be destoyed several times


	virtual fp_t		operator()	(V&  X)			{
		fp_t   y = (*this->wrapped_loft_v)(X); 
								assert(log_file->good());
		// 1 - iter no(gp: splot label);  2 - hight (y);  3 - |X-opt| ignored;  4,5 - X coord;  
		*log_file << format("%d \t %g \t %g \t %22.15g \n")   % (this->wrapped_loft_v->iter())  % this->wrapped_loft_v->opt_distance(X)   %y   %X; 
		return  y;
	};
 };
 }; // namespace lopti
